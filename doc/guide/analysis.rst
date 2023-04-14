.. _analysis:

Tools for analysis
==================

This section explains how to use some tools of the package in order to analyse the data. There are certain tools here which are not available for some DFT code interfaces. Please refer to the DFT package interface converter documentation (see :ref:`conversion`) on how to interface the required DFT outputs into the HDF5 files needed for the tools discussed here. This section will assume that the user has converted the required DFT data.

The following routines require a self energy on the real frequency axis if the user specifies the inputs `with_Sigma` and `with_dc`: 

  * :meth:`density_of_states <dft.sumk_dft_tools.SumkDFTTools.density_of_states>` for the momentum-integrated spectral function including self energy effects and
  * :meth:`spaghettis <dft.sumk_dft_tools.SumkDFTTools.spaghettis>` for the momentum-resolved spectral function (i.e. ARPES)
  * :meth:`spectral_contours <dft.sumk_dft_tools.SumkDFTTools.spectral_contours>` for the k-resolved spectral function on a specific k-mesh (i.e., spectral function on a two dimensional k-mesh)

.. note::
  This package does NOT provide an explicit method to do an **analytic continuation** of
  self energies and Green functions from Matsubara frequencies to the real-frequency axis,
  but a list of options available within the TRIQS framework is given :ref:`here <ac>`.
  Keep in mind that all these methods have to be used very carefully!

Otherwise, without these options, the spectral functions from the inputs of the interfaced DFT code will be used. 

The other routines presented here use the Matsubara self-energy.


Initialisation
--------------

All tools described below are collected in an extension of the :class:`SumkDFT <dft.sumk_dft.SumkDFT>` class and are
loaded by importing the module :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`::

  from triqs_dft_tools.sumk_dft_tools import *

The initialisation of the class is equivalent to that of the :class:`SumkDFT <dft.sumk_dft.SumkDFT>`
class::

  SK = SumkDFTTools(hdf_file = filename + '.h5')

Note that all routines available in :class:`SumkDFT <dft.sumk_dft.SumkDFT>` are also available here.

If required, we have to load and initialise the real-frequency self energy. Most conveniently,
you have your self energy already stored as a real-frequency :class:`BlockGf <triqs.gf.BlockGf>` object
in a hdf5 file::

  with HDFArchive('case.h5', 'r') as ar:
        SigmaReFreq = ar['dmft_output']['Sigma_w']

You may also have your self energy stored in text files. For this case the :ref:`TRIQS <triqslibs:welcome>` library offers
the function :meth:`read_gf_from_txt`, which is able to load the data from text files of one Green function block
into a real-frequency :class:`ReFreqGf <triqs.gf.ReFreqGf>` object. Loading each block separately and
building up a :class:´BlockGf <triqs.gf.BlockGf>´ is done with::

  from triqs.gf.tools import *
  # get block names
  n_list = [n for n,nl in SK.gf_struct_solver[0].iteritems()]
  # load sigma for each block - in this example sigma is composed of 1x1 blocks
  g_blocks = [read_gf_from_txt(block_txtfiles=[['Sigma_'+name+'.dat']], block_name=n) for n in n_list]
  # put the data into a BlockGf object
  SigmaReFreq = BlockGf(name_list=n_list, block_list=g_blocks, make_copies=False)

where:
 
  * `block_txtfiles` is a rank 2 square np.array(str) or list[list[str]] holding the file names of one block and
  * `block_name` is the name of the block.

It is important that each data file has to contain three columns: the real-frequency mesh, the real part and the imaginary part
of the self energy - exactly in this order! The mesh should be the same for all files read in and non-uniform meshes are not supported.

Finally, we set the self energy into the `SK` object::

    SK.set_Sigma([SigmaReFreq])

and additionally set the chemical potential and the double counting correction from the DMFT calculation::

  chemical_potential, dc_imp, dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
  SK.set_mu(chemical_potential)
  SK.set_dc(dc_imp,dc_energ)


Density of states
-----------------

For plotting the density of states, you type::

  SK.density_of_states(mu, broadening, mesh, with_Sigma, with_dc, proj_type, dosocc, save_to_file)

where a brief description of all of the inputs are given in :meth:`density_of_states <dft.sumk_dft_tools.SumkDFTTools.density_of_states>`, which a more in depth discussion of using this routine is given here.

Input parameter breakdown:
  * `mu` (Optional, default is None): The `mu` input specifies the chemical potential to be used in the calculation. By default, this is automatically set to the chemical potential within the SK object.
  * `broadening` (Optional, default 0.01 eV). Traditionally, a small imaginary part is included for a non-interacting spectral function to avoid any numerical artifacts. The user can specify this imaginary part with the `broadening` input. The default value of `0.01 eV` if there is no `broadening` or `mesh` input.
  * `mesh` (Optional, default is None): The `mesh` input requires a real frequency MeshType variable to specificy the frequency on which the spectral function is evaluated on. If the real frequency self-energy is used, then the associated mesh of the self-energy will be used instead. However, this `mesh` input will superseed the SK.mesh real frequency MeshType if it exists.
  * `with_Sigma` and `with_dc` (Optional, both set to True by default): These Boolean options determine whether to include the real frequency self-energy (`with_Sigma`) and double counting term (`with_dc`) from the DMFT calculation. To calculate the DFT+DMFT DOS, both options need to be set to True. If both are set to false, then the DFT DOS will be calculated instead.
  * `proj_type` (Optional, default is None): The type of projection used for the orbital-projected DOS is specified by the optional `proj_type`. These projected spectral functions will be determined alongside the total spectral function. By default, no projected DOS type will be calculated (the corresponding projected arrays will be empty). A full description of the DFT-code dependent options is given below.
  * `dosocc` (Optional, default is False): The occupied part of the spectral function can be generated by setting `dosocc` to True. As a prerequisite to using this option, the user must first calculate the band-resolved density matrix (or occupations) by calling :meth:`occupations <dft.sumk_dft_tools.SumkDFTTools.occupations>` with the Matsubara self-energy.
  * `save_to_file` (Optional, default is True): To save the generated DOS from this routine, then the optional `save_to_file` must be set to True.

Variables returned from this function:
  * DOS - this is a dictionary of numpy arrays with the form of `DOS[spn][n_om]` where `spn` speficies the spin type of the calculation ('`up`, `down`, or combined `ud` which relates to calculations with spin-orbit coupling) and `n_om` is the number of real frequencies as specified by the real frequency MeshType used in the calculation. This array gives the total density of states
  * DOSproj - A dictionary of numpy arrays with the form of `DOSproj[n_shells][spn][n_om]` where `spn` and `n_om` are given above and `n_shells` are the total number of correlated or uncorrelated shells (depending on the input `dos_type`). This array gives the trace of the orbital-projected density of states. 
  * DOSproj_orb - A dictionary of numpy arrays with the form of `DOSproj_orb[n_shells][spn][n_om,dim,dim]` where the variables are defined as given previously with `dim` specifying the orbital dimension of the correlated/uncorrelated shell (depending on the input `dos_type`). This array gives the orbital-projected density of states. 
  

The output files (if `save_to_file = True`) have two (three in the orbital-resolved case) columns representing the frequency and real part of the DOS (and imaginary part of the DOS) in that order. 

The output files are as follows:
  * `DOS_(sp).dat`: The total DOS, where `(dos_type)` and `(sp)` have been defined previously.
  * `DOS_(proj_type)_(sp)_proj(i).dat`: The DOS projected to an orbital with index `(i)`. The index `(i)` refers to tprojindices given in ``SK.shells`` (or ``SK.corr_shells`` for `proj_type` = "wann").
  * `DOS_(proj_type)_(sp)_proj(i)_(m)_(n).dat`: As above, but printed as orbitally-resolved matrix in indices `(m)` and `(n)`. For `d` orbitals, it gives the DOS separately for, e.g., :math:`d_{xy}`, :math:`d_{x^2-y^2}`, and so on.


The different DFT code interfaces read in different variables needed to generate other options for the string type `proj_type` are as follows:  
  * "wann" - calculates the Wannier projected DOS. 
  * "vasp" - calculates the projected DOS from the projectors generated within Vasp.
  * "wien2k" - calculates the projected DOS from the theta projectors generated from dmftproj used to interface Wien2k with TRIQS.

.. image:: images_scripts/DFT_Tools_SVO_DFT_DOS.png
    :width: 600
    :align: center

The figure above shows the DFT SrVO\ :sub:`3`\  density of states generated from 2925 k-points in the irreducible Brillouin zone with the V t\ :sub:`2g`\  Wanner projectors generated within a correlated energy window of [-13.6, 13.6] eV. The `broadening` input has been set to the temperature (i.e., 1/Beta). The total, V t\ :sub:`2g`\  Wannier and occupied total density of states generated from the SK.density_of_states() routine are shown. Note that the noise in the density of states comes from the number of k-points used. This can be removed upon by either using more k-points or using a larger `broadening` value.

Band resolved density matrices
------------------------------

Calculates the band resolved density matrices (occupations) from the Matsubara frequency self-energy.
This is done by calling the following::

  SK.occupations(mu, with_Sigma, with_dc, save_occ):

This is required to generate the occupied DOS in SK.density_of_states() when dosocc is set to True. The `save_occ` optional input (True by default) saves these density matrices to the HDF5 file within the misc_data subgroup. The other variables are the same as defined above. See :meth:`occupations <dft.sumk_dft_tools.SumkDFTTools.occupations>`




Momentum resolved spectral function (with real-frequency self energy)
---------------------------------------------------------------------

Another quantity of interest is the calculated momentum-resolved spectral function A(k, :math:`\omega`) or (correlated) band structure which can directly be compared to ARPES experiments. 
First we have generate the required files from the DFT code of choice and interface them with DFT_Tools, see the guides of the DFT converters (:ref:`conversion`) on how to do this.

This spectral function is calculated by typing::

  SK.spaghettis(mu, broadening, mesh, plot_shift, plot_range, shell_list, with_Sigma, with_dc, proj_type, save_to_file)

Input variables such as `mu`, `broadening` and so on have the same definition as given previously.

The additional parameters are defined as follows:
  * `plot_shift` (Optional, default is 0.0): An additional shift added as `(ik-1)*plot_shift`, where `ik` is the index of the `k` point. This is useful for plotting purposes. 
  * `plot_range` (Optional, default is None): A list with two entries, :math:`\omega_{min}` and :math:`\omega_{max}`, which set the plot range for the output. By default, the full momentum range as given in the self_energy.mesh, input mesh, or SK.mesh is used.
  * `shell_list` (Optional, default is None): A list of integers denoting the shell indices onto which the spectral function is projected if `proj_type` is not None. The resulting function is saved in    the files. By default, the projected spectral function is generated for all shells if `proj_type` is not None. Note for experts: The spectra are not rotated to the local coordinate system used in Wien2k.
  * `proj_type` : This has the same definition as above but with only the "wann" and "wien2k" options implemented.

Variables returned from this function:
  * `Akw` - this is a dictionary of numpy arrays with the form of `Akw[spn][n_ik, n_om]` where `n_k` is the number of k-points used in the calculation. This array gives the total A(k, :math:`\omega`).
  * `pAkw` - A dictionary of numpy arrays with the form of `pAkw[n_shells][spn][n_k, n_om]` where `spn`, `n_om` and `n_shells` have been defined previously. This array gives the trace of the orbital-projected A(k, :math:`\omega`).
  * `pAkw_orb` - A dictionary of numpy arrays with the form of `pAkw_orb[n_shells][spn][n_om,dim,dim]` where the variables have been previously defined. This array gives the orbital-projected A(k, :math:`\omega`).


The output files (if `save_to_file = True`) have three columns representing the k-point index, frequency and A(k, :math:`\omega`) in that order. The output files are as follows:
  * `Akw_(sp).dat`: The total A(k, :math:`\omega`).
  * `Akw_(dos_type)_(sp)_proj(i).dat`: The A(k, :math:`\omega`) projected to shell with index `(i)` (same as defined in the DOS output files). 
  * `Akw_(dos_type)_(sp)_proj(i)_(m)_(n).dat`: As above, but printed as orbitally-resolved matrix in indices `(m)` and `(n)`. 

Also see :meth:`spaghettis <dft.sumk_dft_tools.SumkDFTTools.spaghettis>`.

.. image:: images_scripts/DFT_Tools_SVO_DFT_spaghettis.png
    :width: 1000
    :align: center

The figure above shows the DFT SrVO\ :sub:`3`\  spaghetti plot (generated using V t\ :sub:`2g`\  Wanner projectors generated within a correlated energy window of [-13.6, 13,6] eV). As before, the broadening input has been set to the temperature (i.e., 1/Beta). The left panel shows the total A(k, :math:`\omega`) whereas the right gives the Wannier A(k, :math:`\omega`), both generated from this SK.spaghettis().


Energy contours of the k-resolved Spectral function
---------------------------------------------------

Currently, this has only been implemented for Elk DFT inputs only.

This routine calculates the k-resolved spectral function evaluated at the Fermi level or several energy contours on the k-mesh defined in the converter stage::

  SK.spectral_contours(mu, broadening, mesh, plot_range, FS, with_Sigma, with_dc, proj_type, save_to_file)

Again, certain inputs in this routine have already been defined. Note that only proj_type="wann" has been implemented for this routine. 

Below gives a description of the other inputs used here:
  * `FS` (Optional, default is True): If True, this outputs the spectral function evaluated at the Fermi level (or the frequency value closest to 0.0 eV in the given mesh). If false, then the spectral function is evaluated at the frequency values within the frequency mesh used (or within `plot_range` if specified).

This routine returns the following:
  * `Akw` - A dictionary of numpy arrays with the form of `Akw[spn][n_ik, n_om]` where `spn`, `n_k` and `n_om` being defined as given above. This array gives the total A(k, :math:`\omega`).
  * `pAkw` - A dictionary of numpy arrays with the form of `pAkw[n_shells][spn][n_k, n_om]` (`n_shells` defined above). This array gives the trace of the orbital-projected A(k, :math:`\omega`).
  * `pAkw_orb` - A dictionary of numpy arrays with the form of `pAkw_orb[n_shells][spn][n_om,dim,dim]` (`dim` defined above). This array gives the orbital-projected A(k, :math:`\omega`).


The output files will have the same extensions of `_(sp)_proj(i)_(m)_(n).dat`, `_(sp)_proj(i)_(m)_(n).dat` and `_(sp)_proj(i)_(m)_(n).dat` as described in the Spaghettis routine. 

The files are prepended with either of the following: 
  * For `FS` set to True : The output files' name are prepended with `Akw_FS` and these files contain four columns which are the `k`\ :sub:`x`\, `k`\ :sub:`y`\, `k`\ :sub:`z`\, and `Akw`. Here (`k`\ :sub:`x`\, `k`\ :sub:`y`\, `k`\ :sub:`z`\) are the cartesian reciprocal coordinates, 
  * For `FS` set to False : The output files' name are prepended with `Akw_omega_(iom)` (with `iom` being the frequency mesh index). These files also contain four columns as described above along with a comment at the top of the file which gives the :math:`\omega` frequency value at which the spectral function was evaluated. 

Also see :meth:`spectral_contours <dft.sumk_dft_tools.SumkDFTTools.spectral_contours>`

.. image:: images_scripts/DFT_Tools_SVO_DFT_energy_contours.png
    :width: 1000
    :align: center

The figure above shows the DFT SrVO\ :sub:`3`\  energy contour plots (again, generated using V t\ :sub:`2g`\  Wanner projectors generated within a correlated energy window of [-13.6, 13,6] eV and broadening of 1/Beta). Both panels have been generated on a k-mesh within the first Brilluoin zone on the k\ :sub:`z`\ = 0.0 plane centered at the :math:`\Gamma` point. Here, each panel generated using the outputs from this SK.spectral_contours_plot() routine shows the A(k, :math:`\omega`) evaluated at :math:`\omega` = -0.5 eV (left) and the Fermi level, :math:`\omega` = 0.0 eV, (right).


Partial charges
---------------

Currently, this has only been implemented for Wien2k DFT inputs only.

Since we can calculate the partial charges directly from the Matsubara Green functions, we also do not need a
real-frequency self energy for this purpose. The calculation is done by::

  SK.set_Sigma(SigmaImFreq)
  dm = SK.partial_charges(beta=40.0, with_Sigma=True, with_dc=True)

which calculates the partial charges using the self energy, double counting, and chemical potential as set in the 
`SK` object. On return, `dm` is a list, where the list items correspond to the density matrices of all shells
defined in the list `SK.shells`. This list is constructed by the Wien2k converter routines and stored automatically
in the hdf5 archive. For the structure of `dm`, see also :meth:`partial charges <dft.sumk_dft_tools.SumkDFTTools.partial_charges>`.


