.. _welcome: 

#########
DFTTools
#########

.. sidebar:: DFTTools |PROJECT_VERSION|

   This is the homepage of DFTTools |PROJECT_VERSION|
   For changes see the :ref:`changelog page <changelog>`.
      
      .. image:: _static/logo_github.png
         :width: 75%
         :align: center
         :target: https://github.com/triqs/dft_tools


This :ref:`TRIQS-based <triqslibs:welcome>`-based application is aimed
at ab-initio calculations for 
correlated materials, combining realistic DFT band-structure
calculations with the dynamical mean-field theory. Together with the
necessary tools to perform the DMFT self-consistency loop for
realistic multi-band problems. The package provides a full-fledged
charge self-consistent interface to the `Wien2K package
<http://www.wien2k.at>`_, and `VASP package <https://www.vasp.at>`_. 
In addition, it provides a generic interface for one-shot DFT+DMFT 
calculations, where only the single-particle Hamiltonian in 
orbital space has to be provided. The Hamiltonian can be  
generated from the above mentioned DFT codes, 
`wannier90 <http://www.wannier.org/>`_ output files, or with the 
built-in generic H(k) converter.

Learn how to use this package in the :ref:`documentation` and the :ref:`tutorials`.

.. image:: _static/logo_cea.png
   :width: 14%
   :target: http://ipht.cea.fr

.. image:: _static/logo_x.png
   :width: 14%
   :target: "https://www.cpht.polytechnique.fr

.. image:: _static/logo_cnrs.png
   :width: 14%
   :target: https://www.cnrs.fr

.. image:: _static/logo_erc.jpg
   :width: 14%

.. image:: _static/logo_flatiron.png
   :width: 20%
   :target: https://www.simonsfoundation.org/flatiron

.. image:: _static/logo_simons.jpg
   :width: 20%
   :target: https://www.simonsfoundation.org

    
.. toctree::
   :maxdepth: 2
   :hidden:

   install
   documentation
   tutorials
   issues
   ChangeLog.md
   about
