.. _blockstructure:

Manipulating the Greens functions block structure
=================================================

The DFTTools package includes the general :class:`BlockStructure <dft.block_structure.BlockStructure>` class for manipulating the blocks of Greens functions. In the following, we will introduce its basic and most commonly used functionalities that might show up in an actual DFT+DMFT calculation, and will illustrate them on a very basic fictitious problem.

SumK and Solver structure
-------------------------

The main idea is to have two structures for the Greens functions available. The first one is used in the procedures of the :class:`SumkDFT <dft.sumk_dft.SumkDFT>` to calculate Dysons equations, lattice Greens functions, and so on. The second structure is the one which is used for the solution of the Anderson impurity problem. As a matter of fact, these two structure need not necessarily be the same.

[TODO] Describe the mapping. Discuss why small blocks are important for the QMC solver. Do the example.

Picking orbitals
----------------

In some cases it might happen that for the projection to localised orbitals a full d or f-shells has to be used. However, for the Anderson impurity problem, just a subset of the orbitals are needed. This is the case, e.g., when the projection leads to completely empty or full orbitals that you don't want to include in the AIM. 

For the example we are dealing with it means that .... [CONTINUE]

Basis rotations
---------------

In cases where the Greens function or the local Hamiltonian shows off diagonal entries in the chosen basis, it is often beneficial to rotate to a different basis. This is of particular interest when using a QMC solver, since off-diagonal contributions lead to a famous fermionic sign problem. The :class:`BlockStructure <dft.block_structure.BlockStructure>` class includes methods to perform such basis rotations.

[TODO] Show it on the example

Dismissing off-diagonals
------------------------

As said above, off diagonal contributions lead to some troubles. However,
When you are exactly sure that you know what you are doing, there is functionality to take only the diagonal parts into account in the block structure. Be careful, there is no automatic check whether this approximation is justified or not!

[TODO] Show the example  