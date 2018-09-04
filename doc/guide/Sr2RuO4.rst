.. _Sr2RuO4:

Spin-orbit coupled calculations (single-shot)
=============================================

There are two main ways of including the spin-orbit coupling (SO) term into
DFT+DMFT calculations:

-   by performing a DFT calculation including SO and then doing a DMFT calculation on top, or
-   by performing a DFT calculation without SO and then adding the SO term on the model level.

Treatment of SO in DFT
----------------------

For now, TRIQS/DFTTools does only work with Wien2k when performing calculations with SO.
Of course, the general Hk framework is also possible.
But the way VASP treats SO is fundamentally different to the way Wien2k treats it and the interface does not handle that at the moment.

Therefore, this guide assumes that Wien2k is being used.

First, a Wien2k calculation including SO has to be performed.
For details, we refer the reader to the documentation of Wien2k.
The interface to Wien2k only works when the DFT calculation is done both spin-polarized and with SO (that means that you have to initialize the Wien2k calculation accordingly and then run with ``runsp -sp``).

Performing the projection
~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the final ``x lapw2 -almd -so -up`` and ``x lapw2 -almd -so -dn`` have to be run *on a single core*, which implies that, before, ``x lapw1 -up``, ``x lapw2 -dn``, and ``x lapwso -up`` have to be run in single-core mode (once).

In the ``case.indmftpr`` file, the spin-orbit flag has to be set to ``1`` for the correlated atoms.
For example, for the compound Sr\ :sub:`2`\ RuO\ :sub:`4`, with the struct file :download:`Sr2RuO4.struct <Sr2RuO4/Sr2RuO4.struct>`, we would e.g. use the ``indmftpr`` file :download:`found here <Sr2RuO4/Sr2RuO4.indmftpr>`.
Then, ``dmftproj -sp -so`` has to be called.
As usual, it is important to check for warnings (e.g., about eigenvalues of the overlap matrix) in the output of ``dmftproj`` and adapt the window until these warnings disappear.

Note that in presence of SO, it is not possible to project only onto the :math:`t_{2g}` subshell because it is not an irreducible representation.
A redesign of the orthonormalization procedure might happen in the long term, which might allow that.

We strongly suggest using the :py:meth:`.dos_wannier_basis` functionality of the :py:class:`.SumkDFTTools` class (see :download:`calculate_dos_wannier_basis.py <Sr2RuO4/calculate_dos_wannier_basis.py>`) and compare the Wannier-projected orbitals to the original DFT DOS (they should be more or less equal).
Note that, with SO, there are usually off-diagonal elements of the spectral function, which can also be imaginary.
The imaginary part can be found in the third column of the files ``DOS_wann_...``.
