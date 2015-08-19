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

.. note::
  Add the doc for loading the self energy from a data file. We have to provide this option, because
  in general the user won't has it stored in h5 file!!

Initialisation
--------------

All tools described below are collected in an extension of the :class:`SumkDFT` class and are
loaded by importing the module :class:`SumkDFTTools`::

  from pytriqs.applications.dft.sumk_dft_tools import *

The initialisation of the class is equivalent to that of the :class:`SumkDFT` 
class::

  SK = SumkDFTTools(hdf_file = filename + '.h5')

Note that all routines available in :class:`SumkDFT` are also available here. 

If required, a self energy is load and initialise in the next step. Most conveniently, 
your self energy is already stored as a real frequency :class:`BlockGf` object 
in a hdf5 file::

  ar = HDFArchive('case.h5', 'a')
  SigmaReFreq = ar['dmft_output']['Sigma_w']
  SK.put_Sigma(Sigma_imp = [SigmaReFreq])

Additionally, the chemical potential and the double counting
correction from the DMFT calculation are set, and the archive is closed again::
  
  chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
  SK.set_mu(chemical_potential)
  SK.set_dc(dc_imp,dc_energ)
  del ar


Density of states of the Wannier orbitals
-----------------------------------------

For plotting the density of states of the Wannier orbitals, you type::

  SK.dos_wannier_basis(broadening=0.03, mesh=[om_min, om_max, n_om], with_Sigma=False, with_dc=False, save_to_file=True)

which produces plots between the real frequencies `om_min` and `om_max`, using a mesh of `n_om` points. The parameter 
`broadening` defines an additional Lorentzian broadening, and has the default value of `0.01` eV. To check the Wannier 
density of states after the projection set `with_Sigma` and `with_dc` to `False`.

Partial charges
---------------

Since we can calculate the partial charges directly from the Matsubara Green's functions, we also do not need a
real frequency self energy for this purpose. The calculation is done by::

  SK.put_Sigma(Sigma_imp = SigmaImFreq)
  dm = SK.partial_charges(beta=40.0 with_Sigma=True, with_dc=True)

which calculates the partial charges using the self energy, double counting, and chemical potential as set in the 
`SK` object. On return, `dm` is a list, where the list items correspond to the density matrices of all shells
defined in the list `SK.shells`. This list is constructed by the Wien2k converter routines and stored automatically
in the hdf5 archive. For the detailed structure of `dm`, see the reference manual.

Correlated spectral function (with real frequency self energy)
--------------------------------------------------------------

With this self energy, we can now execute::

  SK.dos_parproj_basis(broadening=broadening)

This produces both the momentum-integrated (total density of states or DOS) and orbitally-resolved (partial/projected DOS) spectral functions.
The variable `broadening` is an additional Lorentzian broadening applied to the resulting spectra.
The output is printed into the files

  * `DOScorr(sp).dat`: The total DOS, where `(sp)` stands for `up`, `down`, or combined `ud`. The latter case
    is relevant for calculations including spin-orbit interaction.
  * `DOScorr(sp)_proj(i).dat`: The DOS projected to an orbital with index `(i)`. The index `(i)` refers to 
    the indices given in ``SK.shells``.
  * `DOScorr(sp)_proj(i)_(m)_(n).dat`: As above, but printed as orbitally-resolved matrix in indices 
    `(m)` and `(n)`. For `d` orbitals, it gives the DOS seperately for, e.g., :math:`d_{xy}`, :math:`d_{x^2-y^2}`, and so on.

Momentum resolved spectral function (with real frequency self energy)
---------------------------------------------------------------------

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
