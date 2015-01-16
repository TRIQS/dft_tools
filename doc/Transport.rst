.. index:: Transport

.. _Transport:

Transport calculations
======================

.. index:: Theory

Formalism
---------
The conductivity and the Seebeck coefficient in direction :math:`\alpha\beta` are defined as [#transp]_:

.. math::

   \sigma_{\alpha\beta} = \beta e^{2} A_{0,\alpha\beta}  \ \ \  \text{and} \ \ \  S_{\alpha\beta} = -\frac{k_B}{|e|}\frac{A_{1,\alpha\beta}}{A_{0,\alpha\beta}}, 

in which the kinetic coefficients :math:`A_{n,\alpha\beta}` are given by

.. math::
  
   A_{n,\alpha\beta} = N_{sp} \pi \hbar \int{d\omega \left(\beta\omega\right)^n f\left(\omega\right)f\left(-\omega\right)\Gamma_{\alpha\beta}\left(\omega\right)}.

Here :math:`N_{sp}` is the spin factor and :math:`f(\omega)` is the Fermi function. The transport distribution :math:`\Gamma_{\alpha\beta}\left(\omega\right)` is defined as

.. math::
  
   \Gamma_{\alpha\beta}\left(\omega\right) = \frac{1}{V} \sum_k Tr\left(v_{k,\alpha}A_{k}(\omega)v_{k,\beta}A_{k}\left(\omega\right)\right),

where :math:`V` is the unit cell volume. In multi-band systems the velocities :math:`v_{k}` and the spectral function :math:`A(k,\omega)` are matrices in the band indices :math:`i` and :math:`j`.
The frequency depended optical conductivity is given by

.. math::

   \sigma(\Omega) = N_{sp} \pi e^2 \hbar \int{d\omega \Gamma_{\alpha\beta}(\omega+\Omega/2,\omega-\Omega/2)\frac{f(\omega-\Omega/2)-f(\omega+\Omega/2)}{\Omega}}.

.. index:: Prerequisite

Prerequisites
-------------
First perform a standard DFT+DMFT calculation for your desired material (see :ref:`DFTDMFTmain`) and obtain the real-frequency self energy by doing an
analytic continuation.

.. note::
   It is crucial to perform the analytic continuation in such a way that the obtained real-frequency self energy is accurate around the Fermi energy as only its
   low energy structure influences the final results!

Besides the self energy the Wien2k files read by the transport converter are:
   * :file:`.struct`: The lattice constants specified in the struct file are used to calculate the unit cell volume.
   * :file:`.outputs`: In this file the k-point symmetries are given.
   * :file:`.oubwin`: Contains the indices of the bands within the projected subspace (written by :program:`dmftproj`) for each k-point.
   * :file:`.pmat`: This file is the output of the Wien2k optics package and contains the velocity (momentum) matrix elements between all bands in the desired energy
     window for each k-point. How to use the optics package is described below.
   * :file:`.h5`: The hdf file has to be present and should contain the dft_input subgroup. Otherwise :class:`convert_dft_input` needs to be called before :class:`convert_transport_input`.

.. index:: Wien2k optics package

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
Here the energy window can be chosen according to the window used for :program:`dmftprj`. However, keep in mind that energies have to be specified in absolute values! Furthermore it is important 
to set line 6 to ON for printing the matrix elements to the :file:`.pmat` file.

.. index:: Using the transport code.

Using the transport code
------------------------

First we have to read the Wien2k files and store the relevant information in the hdf-file::

    from pytriqs.applications.dft.converters.wien2k_converter import *
    from pytriqs.applications.dft.sumk_dft_tools import *

    Converter = Wien2kConverter(filename='case', repacking=True)
    Converter.convert_transport_input()

    SK = SumkDFTTools(hdf_file='case.h5', use_dft_blocks=True)

Additionally we need to read and set the self energy::

    ar = HDFArchive('case_Sigma.h5', 'a')
    SK.put_Sigma(Sigma_imp = [ar['dmft_transp_output']['Sigma_w']])
    del ar

As next step we can calculate the transport distribution :math:`\Gamma_{\alpha\beta}(\omega)`::

    SK.transport_distribution(directions=['xx'], Om_mesh=[0.00, 0.02], energy_window=[-0.3,0.3], 
                                                 with_Sigma=True, broadening=0.0, beta=40)

The parameters are:
    * `directions`: :math:`\alpha` and :math:`\beta` (e.g. xx, yy, xz, ...)
    * `Om_mesh`: :math:`\Omega`-mesh for the optical conductivity. Note that the code repines this mesh to the closest values on the self energy mesh! The new mesh is stored in `Om_meshr`. 
      The Seebeck coefficient is only calculated if :math:`\Omega=0.0` is included.
    * `energy_window`: Limits for the integration over :math:`\omega`. (Due to the Fermi functions the integrand is only of considerable size in a small 
      window around the Fermi energy.)
    * `with_Sigma`: If this parameter is set to False then Sigma is set to 0 (i.e. the DFT band structure :math:`A(k,\omega)` is taken).
    * `broadening`: The numerical broadening should be set to a finite value for with_Sigma = False.

The resulting transport distribution is not automatically saved, but this can be easily achieved with::
    
    SK.save(['Gamma_w','Om_meshr','omega','directions'])
    SK.load(['Gamma_w','Om_meshr','omega','directions'])  

Finally the optical conductivity :math:`\sigma(\Omega)` and the Seebeck coefficient :math:`S` can be obtained with::

    SK.conductivity_and_seebeck(beta=40)
    SK.save(['seebeck','optic_cond']) 

It is strongly advised to check convergence in the number of k-points!

.. index:: References

References
----------

.. [#transp] `V. S. Oudovenko, G. Palsson, K. Haule, G. Kotliar, S. Y. Savrasov, Phys. Rev. B 73, 035120 (2006) <http://link.aps.org/doi/10.1103/PhysRevB.73.0351>`_
.. [#userguide] `P. Blaha, K. Schwarz, G. K. H. Madsen, D. Kvasnicka, J. Luitz, ISBN 3-9501031-1-2 <http://www.wien2k.at/reg_user/textbooks/usersguide.pdf>`_
