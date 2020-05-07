.. _soc:

Spin-orbit coupled calculations (single-shot)
=============================================

There are two main ways of including the spin-orbit coupling (SOC) term into
DFT+DMFT calculations:

-   by performing a DFT calculation including SOC and then doing a DMFT calculation on top, or
-   by performing a DFT calculation without SOC and then adding the SOC term on the model level.

The second variant is a bit more involved and needs quite some expertise, so this guide will cover only the first variant with SOC included in the DFT calculations.

Treatment of SOC in DFT
-----------------------

For now, TRIQS/DFTTools does only work with Wien2k when performing calculations with SO.
The treatment of SOC in the VASP package is fundamentally different to the way Wien2k treats it, and the interface does not handle that at the moment.
Therefore, this guide describes how to do an SOC calculation using the Wien2k DFT package.

First, a Wien2k calculation including SOC has to be performed.
For details, we refer the reader to the `documentation of Wien2k <http://susi.theochem.tuwien.ac.at/reg_user/textbooks/>`_ . As a matter of fact, we need the output for the DFT band structure for both spin directions explicitly. That means that one needs to do a spin-polarised DFT calculation with SOC, but, however, with magnetic moment set to zero. In the Wien2k initialisation procedure, one can choose for the option -nom when ``lstart`` is called. This means that the charge densities are initialised without magnetic splitting. The SOC calculation is then performed in a standard way as described in the Wien2k manual.

Performing the projection
~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the final ``x lapw2 -almd -so -up`` and ``x lapw2 -almd -so -dn`` have to be run *on a single core*, which implies that, before, ``x lapw2 -up``, ``x lapw2 -dn``, and ``x lapwso -up`` have to be run in single-core mode (at least once).

In the ``case.indmftpr`` file, the spin-orbit flag has to be set to ``1`` for the correlated atoms.
For example, for the compound Sr\ :sub:`2`\ MgOsO\ :sub:`6`, with the struct file :download:`Sr2MgOsO6.struct <Sr2MgOsO6/Sr2MgOsO6.struct>`, we would, e.g., use the ``indmftpr`` file :download:`found here <Sr2MgOsO6/Sr2MgOsO6_SOC.indmftpr>`.
Then, ``dmftproj -sp -so`` has to be called.
As usual, it is important to check for warnings (e.g., about eigenvalues of the overlap matrix) in the output of ``dmftproj`` and adapt the window until these warnings disappear.

Note that in presence of SOC, it is not possible to project only onto the :math:`t_{2g}` subshell because it is not an irreducible representation.


We strongly suggest using the :py:meth:`.dos_wannier_basis` functionality of the :py:class:`.SumkDFTTools` class (see :download:`calculate_dos_wannier_basis.py <Sr2RuO4/calculate_dos_wannier_basis.py>`) and compare the Wannier-projected orbitals to the original DFT DOS (they should be more or less equal).
Note that, with SOC, there are usually off-diagonal elements of the spectral function, which can also be imaginary.
The imaginary part can be found in the third column of the files ``DOS_wann_...``.

After the projection, one can proceed with the DMFT calculation. However, two things need to be noted here. First, since the spin is not a good quantum number any more, there are off-diagonal elements in the hybridisation function and the local Hamiltonian coupling the two spin directions. This will eventually lead to a fermonic sign problem when QMC is used as a impurity solver. Second, although the :math:`e_{g}` subshell needs to be included in the projection, it can in many cases be neglected in the solution of the Anderson impurity model, after a transformation to a rotated local basis is done. This basis, diagonalising the local Hamiltonian in the presence of SOC, is often called the numerical j-basis. How rotations are performed is described in :ref:`basisrotation`, and the cutting of the orbitals in :ref:`blockstructure`.

A DMFT calculation including SOC for Sr2MgOsO6 is included in the :ref:`tutorials`.
