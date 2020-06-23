.. _Transport:

Transport calculations
============================

Formalism
---------
The conductivity, the Seebeck coefficient and the electronic contribution to the thermal conductivity in direction :math:`\alpha\beta` are defined as [#transp1]_ [#transp2]_:

.. math::

   \sigma_{\alpha\beta} = \beta e^{2} A_{0,\alpha\beta}
   
.. math::
    
    S_{\alpha\beta} = -\frac{k_B}{|e|}\frac{A_{1,\alpha\beta}}{A_{0,\alpha\beta}}, 

.. math::

    \kappa^{\text{el}}_{\alpha\beta} = k_B \left(A_{2,\alpha\beta} - \frac{A_{1,\alpha\beta}^2}{A_{0,\alpha\beta}}\right),

in which the kinetic coefficients :math:`A_{n,\alpha\beta}` are given by

.. math::
  
   A_{n,\alpha\beta} = N_{sp} \pi \hbar \int{d\omega \left(\beta\omega\right)^n f\left(\omega\right)f\left(-\omega\right)\Gamma_{\alpha\beta}\left(\omega,\omega\right)}.

Here :math:`N_{sp}` is the spin factor and :math:`f(\omega)` is the Fermi function. The transport distribution :math:`\Gamma_{\alpha\beta}\left(\omega_1,\omega_2\right)` is defined as

.. math::
  
   \Gamma_{\alpha\beta}\left(\omega_1,\omega_2\right) = \frac{1}{V} \sum_k Tr\left(v_{k,\alpha}A_{k}(\omega_1)v_{k,\beta}A_{k}\left(\omega_2\right)\right),

where :math:`V` is the unit cell volume. In multi-band systems the velocities :math:`v_{k}` and the spectral function :math:`A(k,\omega)` are matrices in the band indices :math:`i` and :math:`j`.
The frequency depended optical conductivity is given by

.. math::

   \sigma(\Omega) = N_{sp} \pi e^2 \hbar \int{d\omega \Gamma_{\alpha\beta}(\omega+\Omega/2,\omega-\Omega/2)\frac{f(\omega-\Omega/2)-f(\omega+\Omega/2)}{\Omega}}.


Prerequisites
-------------
First perform a standard :ref:`DFT+DMFT calculation <full_charge_selfcons>` for your desired material and obtain the
real-frequency self energy.

.. note::
  If you use a CT-QMC impurity solver you need to perform an **analytic continuation** of
  self energies and Green functions from Matsubara frequencies to the real-frequency axis!
  This packages does NOT provide methods to do this, but a list of options available within the TRIQS framework
  is given :ref:`here <ac>`. Keep in mind that all these methods have to be used very carefully. Especially for optics calculations
  it is crucial to perform the analytic continuation in such a way that the real-frequency self energy
  is accurate around the Fermi energy as low-energy features strongly influence the final results.

Besides the self energy the Wien2k files read by the transport converter (:meth:`convert_transport_input <dft.converters.wien2k.Wien2kConverter.convert_transport_input>`) are:
   * :file:`.struct`: The lattice constants specified in the struct file are used to calculate the unit cell volume.
   * :file:`.outputs`: In this file the k-point symmetries are given.
   * :file:`.oubwin`: Contains the indices of the bands within the projected subspace (written by :program:`dmftproj`) for each k-point.
   * :file:`.pmat`: This file is the output of the Wien2k optics package and contains the velocity (momentum) matrix elements between all bands in the desired energy
     window for each k-point. How to use the optics package is described below.
   * :file:`.h5`: The hdf5 archive has to be present and should contain the dft_input subgroup. Otherwise :meth:`convert_dft_input <dft.converters.wien2k.Wien2kConverter.convert_dft_input>` needs to be called before :meth:`convert_transport_input <dft.converters.wien2k.Wien2kConverter.convert_transport_input>`.


Wien2k optics package
---------------------

The basics steps to calculate the matrix elements of the momentum operator with the Wien2k optics package are:
    1) Perform a standard Wien2k calculation for your material.
    2) Run `x kgen` to generate a dense k-mesh. 
    3) Run `x lapw1`.
    4) For metals change TETRA to 101.0 in :file:`case.in2`.
    5) Run `x lapw2 -fermi`.
    6) Run `x optic`. 

Additionally the input file :file:`case.inop` is required. A detail description on how to setup this file can be found in the Wien2k user guide [#userguide]_ on page 166.
The optics energy window should be chosen according to the window used for :program:`dmftproj`. Note that the current version of the transport code uses only the smaller
of those two windows. However, keep in mind that the optics energy window has to be specified in absolute values and NOT relative to the Fermi energy! 
You can read off the Fermi energy from the :file:`case.scf2` file. Please do not set the optional parameter NBvalMAX in :file:`case.inop`.
Furthermore it is necessary to set line 6 to "ON" and put a "1" in the following line to enable the printing of the matrix elements to :file:`case.pmat`.


Using the transport code
------------------------

First we have to read the Wien2k files and store the relevant information in the hdf5 archive::

    from triqs_dft_tools.converters.wien2k import *
    from triqs_dft_tools.sumk_dft_tools import *

    Converter = Wien2kConverter(filename='case', repacking=True)
    Converter.convert_transport_input()

    SK = SumkDFTTools(hdf_file='case.h5', use_dft_blocks=True)

The converter :meth:`convert_transport_input <dft.converters.wien2k.Wien2kConverter.convert_transport_input>`
reads the required data of the Wien2k output and stores it in the `dft_transp_input` subgroup of your hdf file. 
Additionally we need to read and set the self energy, the chemical potential and the double counting::

    with HDFArchive('case.h5', 'r') as ar:
        SK.set_Sigma([ar['dmft_output']['Sigma_w']])
        chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
        SK.set_mu(chemical_potential)
        SK.set_dc(dc_imp,dc_energ)

As next step we can calculate the transport distribution :math:`\Gamma_{\alpha\beta}(\omega)`::

    SK.transport_distribution(directions=['xx'], Om_mesh=[0.0, 0.1], energy_window=[-0.3,0.3], 
                                                 with_Sigma=True, broadening=0.0, beta=40)

Here the transport distribution is calculated in :math:`xx` direction for the frequencies :math:`\Omega=0.0` and :math:`0.1`. 
To use the previously obtained self energy we set with_Sigma to True and the broadening to :math:`0.0`.
As we also want to calculate the Seebeck coefficient and the thermal conductivity we have to include :math:`\Omega=0.0` in the mesh. 
Note that the current version of the code repines the :math:`\Omega` values to the closest values on the self energy mesh.
For complete description of the input parameters see the :meth:`transport_distribution reference <dft.sumk_dft_tools.SumkDFTTools.transport_distribution>`.

The resulting transport distribution is not automatically saved, but this can be easily achieved with::
    
    SK.save(['Gamma_w','Om_meshr','omega','directions'])

You can retrieve it from the archive by::
      
    SK.Gamma_w, SK.Om_meshr, SK.omega, SK.directions = SK.load(['Gamma_w','Om_meshr','omega','directions'])  

Finally the optical conductivity :math:`\sigma(\Omega)`, the Seebeck coefficient :math:`S` and the thermal conductivity :math:`\kappa^{\text{el}}` can be obtained with::

    SK.conductivity_and_seebeck(beta=40)
    SK.save(['seebeck','optic_cond','kappa']) 

It is strongly advised to check convergence in the number of k-points!


References
----------

.. [#transp1] `V. S. Oudovenko, G. Palsson, K. Haule, G. Kotliar, S. Y. Savrasov, Phys. Rev. B 73, 035120 (2006) <http://link.aps.org/doi/10.1103/PhysRevB.73.0351>`_
.. [#transp2] `J. M. Tomczak, K. Haule, T. Miyake, A. Georges, G. Kotliar, Phys. Rev. B 82, 085104 (2010) <https://link.aps.org/doi/10.1103/PhysRevB.82.085104>`_
.. [#userguide] `P. Blaha, K. Schwarz, G. K. H. Madsen, D. Kvasnicka, J. Luitz, ISBN 3-9501031-1-2 <http://www.wien2k.at/reg_user/textbooks/usersguide.pdf>`_
