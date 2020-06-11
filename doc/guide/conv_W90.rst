.. _convW90:

Wannier90 Converter
===================

Using this converter it is possible to convert the output of
`wannier90 <http://wannier.org>`_
Maximally Localized Wannier Functions (MLWF) and create a HDF5 archive
suitable for one-shot DMFT calculations with the
:class:`SumkDFT <dft.sumk_dft.SumkDFT>` class.

The user must supply two files in order to run the Wannier90 Converter:

#. The file :file:`seedname_hr.dat`, which contains the DFT Hamiltonian
   in the MLWF basis calculated through :program:`wannier90` with ``hr_plot = true``
   (please refer to the :program:`wannier90` documentation).
#. A file named :file:`seedname.inp`, which contains the required
   information about the :math:`\mathbf{k}`-point mesh, the electron density,
   the correlated shell structure, ... (see below).

Here and in the following, the keyword ``seedname`` should always be intended
as a placeholder for the actual prefix chosen by the user when creating the
input for :program:`wannier90`.
Once these two files are available, one can use the converter as follows::

    from triqs_dft_tools.converters import Wannier90Converter
    Converter = Wannier90Converter(seedname='seedname')
    Converter.convert_dft_input()

The converter input :file:`seedname.inp` is a simple text file with
the following format (do not use the text/comments in your input file):

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

The last line of :file:`seedname.inp` is the DFT Fermi energy (in eV), which is subtracted from the onsite terms in the :file:`seedname_hr.dat` file. This is recommended since some functions in DFTTools  implicitly assume a Fermi energy of 0 eV. 

The converter will analyse the matrix elements of the local Hamiltonian
to find the symmetry matrices `rot_mat` needed for the global-to-local
transformation of the basis set for correlated orbitals
(see section :ref:`hdfstructure`).
The matrices are obtained by finding the unitary transformations that diagonalize
:math:`\langle w_i | H_I(\mathbf{R}=0,0,0) | w_j \rangle`, where :math:`I` runs
over the correlated shells and `i,j` belong to the same shell (more details elsewhere...).
If two correlated shells are defined as equivalent in :file:`seedname.inp`,
then the corresponding eigenvalues have to match within a threshold of 10\ :sup:`-5`,
otherwise the converter will produce an error/warning.
If this happens, please carefully check your data in :file:`seedname_hr.dat`.
This method might fail in non-trivial cases (i.e., more than one correlated
shell is present) when there are some degenerate eigenvalues:
so far tests have not shown any issue, but one must be careful in those cases
(the converter will print a warning message).

The current implementation of the Wannier90 Converter has some limitations:

* Since :program:`wannier90` does not make use of symmetries (symmetry-reduction
  of the :math:`\mathbf{k}`-point grid is not possible), the converter always
  sets ``symm_op=0`` (see the :ref:`hdfstructure` section).
* No charge self-consistency possible at the moment.
* Calculations with spin-orbit (``SO=1``) are not supported.
* The spin-polarized case (``SP=1``) is not yet tested.
* The post-processing routines in the module
  :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`
  were not tested with this converter.
* ``proj_mat_all`` are not used, so there are no projectors onto the
  uncorrelated orbitals for now.


