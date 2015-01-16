
.. module:: pytriqs.applications.dft_tools

.. _documentation:

Documentation
=============

A (not so) quick tour
---------------------

.. toctree::
   :maxdepth: 1

   interface
   DFTDMFTmain
   advanced
   analysis
   selfcons
   Ce-HI
   Transport

Projective Wannier functions: the dmftproj package
--------------------------------------------------

In addition to the python-related modules, TRIQS also
provides the Wien2k add-on :program:`dmftproj`. It takes the
information about the wave functions calculated by the `Wien2k package
<http://www.wien2k.at>`_, and constructs projected Wannier functions
that are used as localised orbitals for the DMFT calculation. 

The program :program:`dmftproj` is written in the flavor of the
`Wien2k package <http://www.wien2k.at>`_ without python
support. A detailed description of the usage and options of
:program:`dmftproj`
can be found in :download:`this extensive tutorial <TutorialDmftproj.pdf>`. In
addition, it contains also a description of the Wien2k scripts that
are necessary to do the full charge self-consistent calculations.

Reference manual
----------------

This is the reference manual for the python routines.

.. toctree::
   :maxdepth: 2

   reference/h5structure
   reference/sumk_dft


