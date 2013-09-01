
.. module:: pytriqs.applications.dft_tools

.. _documentation:

Documentation
=============


.. toctree::
   :maxdepth: 1

   interface
   LDADMFTmain
   advanced
   analysis
   selfcons
   Ce-HI

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

