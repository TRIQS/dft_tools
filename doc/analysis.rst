.. _analysis:

Analysing tools
===============

This section explains how to use some tools of the package in order to analyse the data.

.. warning::
  The package does NOT provide an explicit method to do an analytic continuation of the
  self energies and Green functions from Matsubara frequencies to the real frequancy axis! 
  There are methods included e.g. in the ALPS package, which can be used for these purposes. But
  be careful: All these methods have to be used very carefully!!

The analysing tools can be found in an extension of the :class:`SumkLDA` class, they are
loaded by::

  from pytriqs.applications.dft.sumk_lda_tools import *

This import the module ``SumkLDATools``. There are two practical tools, for which you don't
need a self energy on the real axis:

  * The density of states of the Wannier orbitals.
  * Partial charges according to the Wien2k definition.

Other routines need the self energy on the real frequency axis. If you managed to get them, you can
calculate

  * the momentum-integrated spectral function including self-energy effects.
  * the momentum-resolved spectral function (i.e. ARPES)

The initialisation of the class is completely equivalent to the initialisation of the :class:`SumkLDA` 
class::

  SK = SumkLDATools(hdf_file = filename)

By the way, all routines available in :class:`SumkLDA` are also available here. 

Routines without real-frequency self energy
-------------------------------------------

For plotting the 
density of states of the Wannier orbitals, you simply type::

  SK.check_input_dos(om_min, om_max, n_om)

which produces plots between real frequencies `ommin` and `ommax`, using a mesh of `N_om` points. There
is an optional parameter, `broadening`, which defines an additional Lorentzian broadening, and is set to `0.01` 
by default.

Since we can calculate the partial charges directly from the Matsubara Green's functions, we also don't need a
real frequency self energy for this purpose. The calculation is done by::

  ar = HDFArchive(SK.hdf_file)
  SK.put_Sigma([ ar['SigmaF'] ])
  del ar
  dm = SK.partial_charges()

which calculates the partial charges using the input that is stored in the hdf5 file (self energy, double counting, 
chemical potential). Here we assumed that the final self energy is stored as `SigmaF` in the archive. 
On return, dm is a list, where the list items correspond to the density matrices of all shells
defined in the list ``SK.shells``. This list is constructed by the Wien2k converter routines and stored automatically
in the hdf5 archive. For the detailed structure of `dm`, see the reference manual.


Routines with real-frequency self energy
----------------------------------------

In order to plot data including correlation effects on the real axis, one has to provide the real frequency self energy. 
Most conveniently, it is stored as a real frequency :class:`BlockGf` object in the hdf5 file. There is one important thing to
keep in mind. The real frequency self energy has to carry the note `ReFreq`::

  SigmaReFreq.note = 'ReFreq'

This tells the SumkLDA routines, that it is indeed a real frequency Greens function. Supposed you have your self energy now
in the archive, you can type::

  ar=HDFArchive(SK.hdf_file)
  SK.put_Sigma([ ar['SigmaReFreq'] ])
  del ar

This loads the self energy and puts it into the :class:`SumkLDA` class. The chemical potential as well as the double
counting correction was already read in the initialisation process.

With this self energy, we can do now::

  SK.dos_partial()

This produces the momentum-integrated spectral functions (density of states, DOS), also orbitally resolved. 
The output is printed into the files

  * `DOScorr(sp).dat`: The total DOS. `(sp)` stands for `up`, `down`, or combined `ud`. The latter case
    is relevant for calculations including spin-orbit interaction.
  * `DOScorr(sp)_proj(i).dat`: The DOS projected to an orbital with index `(i)`. The index `(i)` refers to 
    the indices given in ``SK.shells``.
  * `DOScorr(sp)_proj(i)_(m)_(n).dat`: Sames as above, but printed as orbitally resolved matrix in indices 
    `(m)` and `(n)`. For `d` orbitals, it gives separately the DOS
    for, e.g., :math:`d_{xy}`, :math:`d_{x^2-y^2}`, and so on.

Another quantity of interest is the momentum-resolved spectral function, which can directly be compared to ARPES
experiments. We assume here that we already converted the output of the :program:`dmftproj` program with the 
converter routines, see :ref:`interfacetowien`. The spectral function is calculated by::

  SK.spaghettis(broadening)

The variable `broadening`1 is an additional Lorentzian broadening that is added to the resulting spectra. The output is
written as the 3-column files ``Akw(sp).dat``, where `(sp)` has the same meaning as above. The output format is 
`k`, :math:`\omega`, `value`. Optional parameters are

  * `shift`: An additional shift, added as `(ik-1)*shift`, where `ik` is the index of the `k` point. Useful for plotting purposes, 
    standard value is 0.0.
  * `plotrange`: A python list with two entries, first being :math:`\omega_{min}`, the second :math:`\omega_{max}`, setting the plot
    range for the output. Standard value is `None`, in this case the momentum range as given in the self energy is plotted. 
  * `ishell`: If this is not `None` (standard value), but an integer, the spectral function projected to the orbital with index `ishell` 
    is plotted to the files. Attention: The spectra are not rotated to the local coordinate system as used in the :program:`Wien2k` 
    program (For experts). 
 

