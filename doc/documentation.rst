.. module:: triqs_dft_tools

.. _documentation:

Documentation
=============

Basic notions
-------------

.. toctree::
   :maxdepth: 1

   basicnotions/first
   basicnotions/dft_dmft
   basicnotions/structure


Construction of local orbitals from DFT
---------------------------------------

.. toctree::
   :maxdepth: 2

   guide/conversion
   h5structure


DFT+DMFT
--------

.. toctree::
   :maxdepth: 2

   guide/dftdmft_singleshot
   guide/dftdmft_selfcons

Advanced Topics
---------------

.. toctree::
   :maxdepth: 1

   guide/blockstructure
   guide/BasisRotation
   guide/soc

Postprocessing
--------------

.. toctree::
   :maxdepth: 1

   guide/analysis
   guide/transport


Reference manual
----------------

This is the reference manual for the python routines.

.. autosummary::
   :recursive:
   :toctree: _python_api
   :template: autosummary_module_template.rst

   block_structure
   converters
   sumk_dft
   sumk_dft_tools
   symmetry
   trans_basis



FAQs
----

.. toctree::
   :maxdepth: 2

   faqs/faqs
