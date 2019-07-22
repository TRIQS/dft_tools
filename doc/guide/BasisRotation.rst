.. _plovasp:

Numerical Treatment of the Sign-Problem: Basis Rotations
=======

When performing calculations with off-diagonal complex hybridisation or local Hamiltonian, one is
often limited by the fermionic sign-problem. However, as the sign is no
physical observable, one can try to improve the calculation by rotating
to another basis.

While the choice of basis to perform the calculation in can be chosen
arbitrarily, two choices which have shown good results are the basis
which diagonalizes the local Hamiltonian or the density matrix of then
system.

The transformation matrix can be stored in the :class:`BlockStructure` and the
transformation is automatically performed when using :class:`SumkDFT`'s :meth:`extract_G_loc`
and :meth:`put_Sigma` (see below).


Finding the Transformation Matrix
-----------------

The :class:`TransBasis` class offers a simple mehod to calculate the transformation
matrices to a basis where either the local Hamiltonian or the density matrix
is diagonal::

    from triqs_dft_tools.trans_basis import TransBasis
    TB = TransBasis(SK)
    TB.calculate_diagonalisation_matrix(prop_to_be_diagonal='eal', calc_in_solver_blocks = True)

    SK.block_structure.transformation = [{'ud':TB.w}]



Transforming Green's functions manually
-----------------

One can transform Green's functions manually using the :meth:`convert_gf` method::

    # Rotate a Green's function from solver-space to sumk-space
    new_gf = block_structure.convert_gf(old_gf, space_from='solver', space_to='sumk')



Automatic transformation during the DMFT loop
-----------------

During a DMFT loop one is switching back and forth between Sumk-Space and Solver-Space
in each iteration. Once the block_structure.transformation property is set, this can be
done automatically::

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
  One must not forget to also transform the interaction Hamiltonian to the diagonal basis!
  This can be done with the :meth:`transform_U_matrix` method. However, due to different 
  conventions in this method, one must pass the conjugated version of the transformation matrix::
  
  U_trans = transform_U_matrix(U, SK.block_structure.transformation[0]['ud'].conjugate())
