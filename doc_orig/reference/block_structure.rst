Block Structure
===============

The `BlockStructure` class allows to change and manipulate
Green functions structures and mappings from sumk to solver.

The block structure can also be written to and read from HDF files.

.. warning::

    Do not write the individual elements of this class to a HDF file,
    as they belong together and changing one without the other can
    result in unexpected results. Always write the BlockStructure
    object as a whole.

    Writing the sumk_to_solver and solver_to_sumk elements
    individually is not implemented.

.. autoclass:: triqs_dft_tools.block_structure.BlockStructure
   :members:
   :show-inheritance:
