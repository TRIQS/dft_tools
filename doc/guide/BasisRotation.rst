.. _basisrotation:

Automatic basis rotations in DFT+DMFT
=====================================

When performing calculations with off-diagonal terms in the hybridisation function or in the local Hamiltonian, one is
often limited by the fermionic sign-problem slowing down the QMC calculations significantly. This can occur for instance if the crystal shows locally rotated or distorted cages, or when spin-orbit coupling is included. The examples for this are included in the :ref:`tutorials` of this documentation.

However, as the fermonic sign in the QMC calculation is no
physical observable, one can try to improve the calculation by rotating
to another basis. While this basis can in principle be chosen arbitrarily, 
two choices which have shown good results; name the basis sets that diagonalize the local Hamiltonian or the local density matrix of the
system.

As outlined in section :ref:`blockstructure`, the :class:`BlockStructure` includes all necessary functionalities. While it is possible to manually transform each Green's functions and self energies between the *sumk* and the *solver* basis, this leads to cumbersum code and is discouraged. Instead, in order to facilitate the block-structure manipulations for an actual DFT+DMFT calculation, some of the necessary steps are automatically included automatically. As soon as the 
transformation matrix is stored in the :class:`BlockStructure`, the
transformation is automatically performed when using :class:`SumkDFT`'s :meth:`extract_G_loc`,
:meth:`put_Sigma`, and :meth:`calc_dc` (see below).

Setting up the initial solver structure from DFT
------------------------------------------------

Before the actual calculation one has to specify the *solver* basis structure, in which the impurity problem will be tackled. The different available approximation were introduced in section :ref:`blockstructure`. An important feature of DFTTools is that there is an automatic way to determine the entries of the Green's function matrix that are zero by symmetry, when initialising the class::

    from triqs_dft_tools.sumk_dft import *
    SK = SumkDFT(hdf_file,use_dft_blocks='True')

The flag *use_dft_blocks=True* analysis the local density matrix, determines the zero entries, and sets up a minimal *solver* structure. Alternatively, this step can be achieved by (see the reference manual for options)::

    SK.analyse_block_structure()


Finding the transformation matrix
---------------------------------

The SumkDFT class offers a method that can determine transformation matrices to certain new basis. It is called as follows::

    SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal')

Possible option for *prop_to_be_diagonal* are *eal* (diagonalises the local hamiltonian) or *dm* (diagonalises the local density matrix). This routine stores the transformation matrix in the :class:`SK.block_structure` class, such that it can be used to rotate the basis. 



Automatic transformation during the DMFT loop
---------------------------------------------

During a DMFT loop one is often switching back and forth between the unrotated basis (Sumk-Space) and the rotated basis that is used by the QMC Solver.
Once the SK.block_structure.transformation property is set as shown above, this is
done automatically, meaning that :class:`SumkDFT`'s :meth:`extract_G_loc`, :meth:`put_Sigma`, and :meth:`calc_dc` are doing the transformations by default::

    for it in range(iteration_offset, iteration_offset + n_iterations):
        # every GF is in solver space here
        S.G0_iw << inverse(S.Sigma_iw + inverse(S.G_iw))

        # solve the impurity in solver space -> hopefully better sign
        S.solve(h_int = H, **p)

        # calc_dc(..., transform = True) by default
        SK.calc_dc(S.G_iw.density(), U_interact=U, J_hund=J, orb=0, use_dc_formula=DC_type)
        
        # put_Sigma(..., transform_to_sumk_blocks = True) by default
        SK.put_Sigma([S.Sigma_iw])
        
        SK.calc_mu()

        # extract_G_loc(..., transform_to_solver_blocks = True) by default
        S.G_iw << SK.extract_G_loc()[0]

.. warning::
  Before doing the DMFT self-consistency loop, one must not forget to also transform the 
  interaction Hamiltonian to the diagonal basis!
  This can be also be done with a method of the :class:`BlockStructure` class, 
  namely the :meth:`convert_operator` method. Having set up a Hamiltonian in the *sumk* structure, it can easily 
  be transformed to the *solver* structure (including rotations of basis, picking of orbitals, 
  making matrices diagonal, etc) by::

    H_solver = SK.block_structure.convert_operator(H_sumk)

  We refer to the tutorials on how to set up the Hamiltonian H_sumk in selected cases.
  Note that this transformation might generally lead to complex values in the 
  interaction Hamiltonian. Unless you know better and can make everything real, 
  you should take care of using the correct version of the TRIQS CTQMC solver.