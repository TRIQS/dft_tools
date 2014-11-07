
H5 input file and Converters
============================

All the data is stored using the hdf5 standard, as described also in the documentation of the TRIQS package itself. In order to do a DMFT calculation, using input from DFT applications, a converter is needed on order to provide the necessary data in the hdf5 format. 



hdf5 data format
----------------

In order to be used with the DMFT routines, the following data needs to be provided in the hdf5 file. It contains a lot of information in order to perform DMFT calculations for all kinds of situations, e.g. d-p Hamiltonians, more than one correlated atomic shell, or using symmetry operations for the k-summation. We store all data in subgroups of the hdf5 arxive:

:program:`Main data`: There needs to be one subgroup for the main data of the calculation. The default name of this group is `lda_input`. Its contents are

  * `energy_unit`, numpy.float. The unit of energy used for the calculation

  * `n_k`, numpy.int. The number of k-points used for the BZ integration.

  * `k_dep_projection`, numpy.int. If the dimension of the projection operators depend on the k-point, set to 1. Otherwise set to 0.

  * `SP`, numpy.int. 0 for paramagnetic hamiltonian, 1 for spin polarised hamiltonian.

  * `SO`, numpy.int. 1 if spin-orbit interaction is included, 0 otherwise.

  * `charge_below`, numpy.float. The number of electrons in the crystal below the correlated orbitals. 

  * `density_required`, numpy.float. Required total electron density, important for the determination of the chemical potential. The density of the correlated electrons is then `density_required`-`charge_below`. 

  * `symm_op`, numpy.int. 1 if symmetry operations are used for the BZ sums, 0 if all k-points are directly included in the input.

  * `n_shells`, numpy.int. Number of atomic shells for which post processing is possible. This is `not` the number of correlated orbitals! If you have two correlated atoms in the unit cell that are equivalent, `n_shells` is 2!
 
  * `shells`, double list of numpy.int. First dimension: `n_shells`, second dimension: 4. Information about the atomic shells. For each shell, we give a list of 4 numbers: [index, sort, l, dim]. index is the atom index, sort defines the equivalency of the atoms. For instance, with two equivalent atoms in the unit cell, index runs from 0 to 1, but sort can take only one value 0. l is the angular quantum number, and dim the dimension of the atomic shell.

  * `n_corr_shells`, numpy.int. Number of correlated atomic shells, for which correlations are included. This includes atoms which are equivalent by symmetry. If you have two correlated atoms in the unit cell that are equivalent, `n_corr_shells` is 2! 

  * `corr_shells`, double list of numpy.int. First dimension: `n_corr_shells`, second dimension: 6. Information about the correlated orbitals. For each correlated shell, we give a list of 6 numbers: [index, sort, l, dim, SO, irep]. Similar as for `shells`, index is the atom index, sort defines the equivalency of the atoms. For instance, with two equivalent atoms in the unit cell, index runs from 0 to 1, but sort can take only one value 0. l is the angular quantum number, and dim the dimension of the atomic shell. If spin-orbit is included in the calculation, SO is 1, otherwise 0. irep is a dummy integer, set to 0.

  * `use_rotations`, numpy.int. If local and global coordinate systems are used, this falg is set to 1. Otherwise set to 0.

  * `rot_mat`, list of numpy.array.complex. If `use_rotations` is set, then this contains a list of the rotation matrices. You have to give a rotation matrix for each correlated atomic shell, i.e. the length of the list is `n_corr_shells`, the dimension of the matrices is given by the dimension of the atomic shell (c.f. `corr_shells`). You have to give at least the unity matrix, if no rotations are used!

  * `rot_mat_time_inv`, list of numpy.int, length `n_corr_shells`. This is needed only if `SP` is 1. For each correlated shell, set 1 if the coordinate transformation contains inversion, 0 otherwise. If no rotations are used, or `SP` is 0, give a list of zeros. 

  * `n_reps`, numpy.int. Number of irreducible representation of the correlated shell. For instance, if you plan to use the t2g/eg splitting, then set it to 2.

  * `dim_reps`, list of numpy.int, length `n_reps`. Dimension of the representations. In above example, it is [2,3] for eg and t2g sub sets. 

  * `T`, list of numpy.array.complex. For each `inequivalent` correlated shell, this is a general transformation matrix from the spherical harmonics basis to the basis used for the impurity problem, which is most of the time the real cubic harmonics basis. This matrix is used to calculate the 4-index U matrix.

  * `n_orbitals`, numpy.array.int. Dimension [n_k, SP + 1 - SO]. It contains the number of Bloch bands that are included in the projection window for each k point. The second dimension is 1, unless doing spin-polarised DFT calculations without spin orbit coupling, where it is 2. In that case, the number of included bands might depend on the spin projection up/down.

  * `proj_mat`, numpy.array.complex. Dimension [n_k, SP + 1 - SO, n_corr_shells, Max(corr_shell dimensions), Max(n_orbitals)]. This array contains the projection matrices from Bloch bands to Wannier orbitals, where the last two indices are the matrix indices. For efficient storage reasons, all the matrices have to be of the same size. For instance, if the number of Bloch bands has values 4, 5 and 6, all matrices are initialised with size 6. For the k-points with less bands, only the first entries are used, the rest is just zero.

  * `bz_weights`, numpy.array.float. A one-dimensional array of length n_k containing the weights of the k-points for the k summation.

  * `hopping`, numpy.array.complex. Dimension [n_k, SP + 1 - SO, Max(n_orbitals), Max(n_orbitals)]. This array contains the non-interacting Hamiltonian in matrix form at each k point. Again, similar to the projection matrices, all matrices have to be of the same size. 



:program:`Symmetry operations`: In this subgroup we store all the data for applying the symmetry 
    operations in the DMFT loop (in case you want to use symmetry operations). The default name of this subgroup is `SymmCorr`. This information is needed only if symmetry operations are used to do the k summation. To be continued...


Wien2k Converter
----------------

The dft_tools package comes with a converter to use `Wien2k <http://www.wien2k.at>`_ band structure calculations as input for the DMFT part of the calculation, through the construction of projective Wannier functions. The first step is done by the program :program:`dmftproj`, producing text output files. In the second step, this ouput is read and converted into the hdf5 format, using the python module :class:`Wien2kConverter`.

HERE COMES A LISTING OF THE FUNCTIONS.

General and simple H(k) Converter
---------------------------------

The above described converter of the Wien2k input is quite involved, since Wien2k provides a lot of information, e.g. about symmetry operations, that can be used in the calculation. However, sometimes we want to use a light implementation where the input consists basically only of the Hamiltonian matrix in Wannier basis, given at a grid of k points in the first Brillouin zone. For this purpose, a simple converter is included in the package, called :class:`HkConverter`, which is implemented for the simplest case of paramagnetic DFT calculations without spin-orbit coupling. It reads a simple, easy to construct text file, and produces an archive that can be used for the DMFT calculations. An example input file for a structure with one correlated site with 3 t2g orbitals in the unit cell contains the following:

  10               <- n_k

  1.0              <- density_required

  1                <- n_shells

  1 1 2 3          <- shells, as above: iatom, isort, l, dim

  1                <- n_corr_shells

  1 1 2 3 0 0      <- corr_shells, as above: iatom, isort, l, dim, SO, dummy

  2 2 3            <- n_reps, dim_reps (length 2, because eg/t2g splitting)

After this header, we give the Hamiltonian matrices for al the k-points. for each k-point we give first the matrix of the real part, then the matrix of the imaginary part. The projection matrices are set automatically to unity matrices, no rotations, no symmetry operations are used. That means that the symmetry sub group in the hdf5 archive needs not be set, since it is not used. It is furthermore assumed that all k-points have equal weight in the k-sum. Note that the input file should contain only the numbers, not the comments given in above example.

The Hamiltonian matrices can be taken, e.g., from Wannier90, which contructs the Hamiltonian in a maximally localised Wannier basis.

Note that with this simplified converter, no full charge self consistent calculations are possible!



  





