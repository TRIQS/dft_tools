.. _analysis:

Tools for analysis
==================

This section explains how to use some tools of the package in order to analyse the data.

There are two practical tools for which a self energy on the real axis is not needed, namely:

  * :meth:`dos_wannier_basis` for the density of states of the Wannier orbitals and
  * :meth:`partial_charges` for the partial charges according to the Wien2k definition.

However, a real frequency self energy has to be provided by the user to use the methods:

  * :meth:`dos_parproj_basis` for the momentum-integrated spectral function including self energy effects and
  * :meth:`spaghettis` for the momentum-resolved spectral function (i.e. ARPES)

.. warning::
  This package does NOT provide an explicit method to do an **analytic continuation** of the
  self energies and Green functions from Matsubara frequencies to the real frequency axis! 
  There are methods included e.g. in the :program:`ALPS` package, which can be used for these purposes. 
  Keep in mind that all these methods have to be used very carefully!

Most conveniently, you have your self energy already stored as a real frequency :class:`BlockGf` object 
in a hdf5 file, which can be easily loaded::

  ar = HDFArchive(filename+'.h5','r')
  SigmaReFreq = ar['SigmaReFreq']
  del ar

.. note::
    What happens if one has a self energy only in text files...?

You may also store it in text files. If all blocks of your self energy are of dimension 1x1, you store them in `fname_(block)0.dat` files. Here `(block)` is a block name (`up`, `down`, or combined `ud`). In the case when you have matrix blocks, you store them in `(i)_(j).dat` files, where `(i)` and `(j)` are the orbital indices, in the `fname_(block)` directory.

This self energy is loaded and put into the :class:`SumkDFT` class by the function:: 

  SK.constr_Sigma_real_axis(filename, hdf=True, hdf_dataset='SigmaReFreq',n_om=0)

where:
 
  * `filename`: the name of the hdf5 archive file or the `fname` pattern in text files names as described above,  
  * `hdf`: if `True`, the real-axis self energy will be read from the hdf5 file, otherwise from the text files,
  * `hdf_dataset`: the name of dataset where the self energy is stored in the hdf5 file,
  * `n_om`: the number of points in the real-axis mesh (used only if `hdf=False`).
  
The chemical potential as well as the double counting correction were already read in the initialisation process.

Initialisation
--------------

All tools described below are collected in an extension of the :class:`SumkDFT` class and are
loaded by importing the module :class:`SumkDFTTools`::

  from pytriqs.applications.dft.sumk_dft_tools import *

The initialisation of the class is equivalent to that of the :class:`SumkDFT` 
class::

  SK = SumkDFTTools(hdf_file = filename + '.h5')

Note that all routines available in :class:`SumkDFT` are also available here. 

If required, the real frequency self energy is set with::
  
    SK.put_Sigma(Sigma_imp = [ SigmaReFreq ])

Density of states of the Wannier orbitals
-----------------------------------------

For plotting the 
density of states of the Wannier orbitals, you simply type::

  SK.check_input_dos(om_min, om_max, n_om)

which produces plots between the real frequencies `om_min` and `om_max`, using a mesh of `n_om` points. There
is an optional parameter `broadening` which defines an additional Lorentzian broadening, and has the default value of
`0.01` by default.

Partial charges
---------------

Since we can calculate the partial charges directly from the Matsubara Green's functions, we also do not need a
real-frequency self energy for this purpose. The calculation is done by::

  ar = HDFArchive(SK.hdf_file)
  SK.put_Sigma([ ar['SigmaImFreq'] ])
  del ar
  dm = SK.partial_charges()

which calculates the partial charges using the data stored in the hdf5 file, namely the self energy, double counting, and
chemical potential. Here we assumed that the final self energy is stored as `SigmaImFreq` in the archive. 
On return, `dm` is a list, where the list items correspond to the density matrices of all shells
defined in the list `SK.shells`. This list is constructed by the Wien2k converter routines and stored automatically
in the hdf5 archive. For the detailed structure of `dm`, see the reference manual.

Correlated spectral function (with self energy)
-----------------------------------------------

With this self energy, we can now execute::

  SK.dos_partial(broadening=broadening)

This produces both the momentum-integrated (total density of states or DOS) and orbitally-resolved (partial/projected DOS) spectral functions.
The variable `broadening` is an additional Lorentzian broadening applied to the resulting spectra.
The output is printed into the files

  * `DOScorr(sp).dat`: The total DOS, where `(sp)` stands for `up`, `down`, or combined `ud`. The latter case
    is relevant for calculations including spin-orbit interaction.
  * `DOScorr(sp)_proj(i).dat`: The DOS projected to an orbital with index `(i)`. The index `(i)` refers to 
    the indices given in ``SK.shells``.
  * `DOScorr(sp)_proj(i)_(m)_(n).dat`: As above, but printed as orbitally-resolved matrix in indices 
    `(m)` and `(n)`. For `d` orbitals, it gives the DOS seperately for, e.g., :math:`d_{xy}`, :math:`d_{x^2-y^2}`, and so on.

Momentum resolved spectral function (with self energy)
------------------------------------------------------

Another quantity of interest is the momentum-resolved spectral function, which can directly be compared to ARPES
experiments. We assume here that we already converted the output of the :program:`dmftproj` program with the 
converter routines (see :ref:`conversion`). The spectral function is calculated by::

  SK.spaghettis(broadening)

Optional parameters are

  * `shift`: An additional shift added as `(ik-1)*shift`, where `ik` is the index of the `k` point. This is useful for plotting purposes. 
    The default value is 0.0.
  * `plotrange`: A list with two entries, :math:`\omega_{min}` and :math:`\omega_{max}`, which set the plot
    range for the output. The default value is `None`, in which case the full momentum range as given in the self energy is used. 
  * `ishell`: An integer denoting the orbital index `ishell` onto which the spectral function is projected. The resulting function is saved in 
    the files. The default value is `None`. Note for experts: The spectra are not rotated to the local coordinate system used in :program:`Wien2k`.

The output is written as the 3-column files ``Akw(sp).dat``, where `(sp)` is defined as above. The output format is 
`k`, :math:`\omega`, `value`. 
