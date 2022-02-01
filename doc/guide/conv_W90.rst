.. _convW90:

Interface with Wannier90
========================

This interface allows to convert the output of `wannier90 <http://wannier.org>`_
Maximally Localized Wannier Functions (MLWF) and create a HDF5 archive suitable for DMFT calculations with the
:class:`SumkDFT <dft.sumk_dft.SumkDFT>` class. The tasks are parallelized with MPI.

The converter can be run in two different modes, which are specified with the keyword ``bloch_basis`` in the call::

    from triqs_dft_tools.converters import Wannier90Converter
    Converter = Wannier90Converter(seedname='seedname', bloch_basis=False, rot_mat_type='hloc_diag', add_lambda=None)

Here and in the following, the keyword ``seedname`` should always be intended
as a placeholder for the actual prefix chosen by the user when creating the
input for :program:`wannier90`.

Orbital mode
---------------

In the default mode (``bloch_basis = False``), the Converter writes the Hamiltonian in orbital basis, in which case
the projector functions are trivial identity matrices. The user must supply two files:

#. The file :file:`seedname_hr.dat`, which contains the DFT Hamiltonian
   in the MLWF basis calculated through :program:`wannier90` with ``write_hr = true``
   (please refer to the :program:`wannier90` documentation).
#. A file named :file:`seedname.inp`, which contains the required
   information about the :math:`\mathbf{k}`-point mesh, the electron density,
   the correlated shell structure, ... (see below).

Once these two files are available, one can use the converter as follows::

    Converter.convert_dft_input()

The converter input :file:`seedname.inp` is a simple text file with
the following format:

.. literalinclude:: images_scripts/LaVO3_w90.inp

The example shows the input for the perovskite crystal of LaVO\ :sub:`3`
in the room-temperature `Pnma` symmetry. The unit cell contains four
symmetry-equivalent correlated sites (the V atoms) and the total number
of electrons per unit cell is 8 (see second line).
The first line specifies how to generate the :math:`\mathbf{k}`-point
mesh that will be used to obtain :math:`H(\mathbf{k})`
by Fourier transforming :math:`H(\mathbf{R})`.
Currently implemented options are:

* :math:`\Gamma`-centered uniform grid with dimensions
  :math:`n_{k_x} \times n_{k_y} \times n_{k_z}`;
  specify ``0`` followed by the three grid dimensions,
  like in the example above
* :math:`\Gamma`-centered uniform grid with dimensions
  automatically determined by the converter (from the number of
  :math:`\mathbf{R}` vectors found in :file:`seedname_hr.dat`);
  just specify ``-1``

Inside :file:`seedname.inp`, it is crucial to correctly specify the
correlated shell structure, which depends on the contents of the
:program:`wannier90` output :file:`seedname_hr.dat` and on the order
of the MLWFs contained in it. In this example we have four lines for the
four V atoms. The MLWFs were constructed for the t\ :sub:`2g` subspace, and thus
we set ``l`` to 2 and ``dim`` to 3 for all V atoms. Further the spin-orbit coupling (``SO``)
is set to 0 and ``irep`` to 0.
As in this example all 4 V atoms are equivalent we set ``sort`` to 0. We note
that, e.g., for a magnetic DMFT calculation the correlated atoms can be made
inequivalent at this point by using different values for ``sort``.

The number of MLWFs must be equal to, or greater than the total number
of correlated orbitals (i.e., the sum of all ``dim`` in :file:`seedname.inp`).
If the converter finds fewer MLWFs inside :file:`seedname_hr.dat`, then it
stops with an error; if it finds more MLWFs, then it assumes that the
additional MLWFs correspond to uncorrelated orbitals (e.g., the O-\ `2p` shells).
When reading the hoppings :math:`\langle w_i | H(\mathbf{R}) | w_j \rangle`
(where :math:`w_i` is the :math:`i`-th MLWF), the converter also assumes that
the first indices correspond to the correlated shells (in our example,
the V-t\ :sub:`2g` shells). Therefore, the MLWFs corresponding to the
uncorrelated shells (if present) must be listed **after** those of the
correlated shells.
With the :program:`wannier90` code, this can be achieved by listing the
projections for the uncorrelated shells after those for the correlated shells.
In our `Pnma`-LaVO\ :sub:`3` example, for instance, we could use::

    Begin Projections
     V:l=2,mr=2,3,5:z=0,0,1:x=-1,1,0
     O:l=1:mr=1,2,3:z=0,0,1:x=-1,1,0
    End Projections

where the ``x=-1,1,0`` option indicates that the V--O bonds in the octahedra are
rotated by (approximatively) 45 degrees with respect to the axes of the `Pbnm` cell.

The last line of :file:`seedname.inp` is the DFT Fermi energy (in eV), which is subtracted from the onsite
terms in the :file:`seedname_hr.dat` file. This is recommended since some functions in DFTTools implicitly
assume a Fermi energy of 0 eV.

In the orbital mode the Converter supports the addition of a local spin-orbit term, if the Wannier Hamiltonian
describes a t\ :sub:`2g` manifold. Currently, the correct interaction term is only implemented if the default
orbital order of :program:`wannier90` is maintained, i.e. it is assumed to be
:math:`d_{xz,\uparrow}, d_{yz,\uparrow}, d_{xy,\uparrow}, d_{xz,\downarrow}, d_{yz,\downarrow}, d_{xy,\downarrow}`.
The coupling strength can be specified as ``add_lambda = [lambda_x, lambda_y, lambda_z]``,
representative of the orbital coupling terms perpendicular to :math:`[x, y, z]` i.e. :math:`[d_{yz}, d_{xz}, d_{xy}]`,
respectively. Note that it is required to have ``SO=0`` and ``SP=1``.

For spin-orbit coupled systems (from DFT or with the add_lambda parameter),
the orbitals are resorted internally by the converter to the triqs order so, e.g.,
the order from above would become
:math:`d_{xz,\uparrow}, d_{xz,\downarrow}, d_{yz,\uparrow}, d_{yz,\downarrow}, d_{xy,\uparrow}, d_{xy,\downarrow}`.

Band mode
----------------

If ``bloch_basis = True``, the Converter writes the Hamiltonian in the Kohn-Sham basis that was used to construct
the Wannier functions. The projector functions are then given by the transformation from Kohn-Sham to orbital basis.
Note that to do so :program:`wannier90` must be run with ``write_u_matrices = true``. Additionally to the files
described above, the Converter will require the following files:

#. :file:`seedname_u.mat` (and :file:`seedname_u_dis.mat` if disentanglement was used to construct the Wannier functions.) is read to construct the projector functions.
#. :file:`seedname.eig` is read to get the Kohn-Sham band eigenvalues
#. :file:`seedname.nnkp` is read to obtain the band indices of the orbitals selected for the Wannier Hamiltonian
#. :file:`seedname.wout` is read to get the outer energy window to ensure the correct mapping of the disentanglement

Note that in case of disentanglement the user must set the outer energy window (``dis_win_min`` and ``dis_win_max``) explicitly in :program:`wannier90` with an energy
separation of at least :math:`10^{-4}` to the band energies. This means in particular that one should not use the default energy window to avoid subtle bugs.

Additionally, to keep the dimension of the lattice Green's function reasonable, it is recommendable to use the exclude_bands tag for bands completely outside of the energy window.

The Converter currently works with Quantum Espresso and VASP. Additional files are required for each case to obtain
the Fermi weights:

#. :file:`seedname.nscf.out` for Quantum Espresso (the NSCF run must contain the flag ``verbosity = 'high'``)
#. :file:`OUTCAR` and :file:`LOCPROJ` for VASP

Note that in the band mode the user input of the :math:`k`-mesh and the Fermi energy in :file:`seedname.inp` are ignored, since both quantities
are automatically read from the :program:`wannier90` and DFT output. However, the :math:`k`-mesh parameter still has to be specified to comply with the file format.

Rotation matrix
------------------

The converter will analyse the matrix elements of the local Hamiltonian
to find the symmetry matrices `rot_mat` needed for the global-to-local
transformation of the basis set for correlated orbitals
(see section :ref:`hdfstructure`).
If ``rot_mat_type='hloc_diag'``, the matrices are obtained by finding the unitary transformations that diagonalize
:math:`\langle w_i | H_I(\mathbf{R}=0,0,0) | w_j \rangle`, where :math:`I` runs
over the correlated shells and `i,j` belong to the same shell (more details elsewhere...).
If ``rot_mat_type='wannier'``, the matrix for the first correlated shell per impurity will be identity, defining the reference frame,
while the rotation matrices of all other equivalent shells contain the correct mapping into this reference frame.
If two correlated shells are defined as equivalent in :file:`seedname.inp`,
then the corresponding eigenvalues have to match within a threshold of 10\ :sup:`-5`,
otherwise the converter will produce an error/warning.
If this happens, please carefully check your data in :file:`seedname_hr.dat`.
This method might fail in non-trivial cases (i.e., more than one correlated
shell is present) when there are some degenerate eigenvalues:
so far tests have not shown any issue, but one must be careful in those cases
(the converter will print a warning message and turns off the use of rotation matrices,
which leads to an incorrect mapping between equivalent correlated shells).

Note that in the case of spin-orbit coupling, these rotation matrices in general
mix spin and orbital components.

Current limitations
----------------------------------------------

The current implementation of the Wannier90 Converter has some limitations:

* Since :program:`wannier90` does not make use of symmetries (symmetry-reduction
  of the :math:`\mathbf{k}`-point grid is not possible), the converter always
  sets ``symm_op=0`` (see the :ref:`hdfstructure` section).
* The spin-polarized case (``SP=1``) is neither completely implemented nor tested.
* ``proj_mat_all`` are not used, so there are no projectors onto the
  uncorrelated orbitals for now.
