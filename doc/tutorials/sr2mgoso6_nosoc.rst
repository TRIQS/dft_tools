.. _Sr2MgOsO6_noSOC:

Here we will discuss a calculation where off-diagonal matrix elements show up, and will discuss step-by-step how this calculation can be set up.

The full script for this calculation is also provided here (:download:`Sr2MgOsO6_noSOC.py <images_scripts/Sr2MgOsO6_noSOC.py>`).

Note that we do not include spin-orbit coupling here for pedagogical reasons. For the real material it is necessary to include also SOC.


DFT (Wien2k) and Wannier orbitals
=================================

DFT setup
---------

First, we do a DFT calculation, using the Wien2k package. As main input file we have to provide the so-called struct file :file:`Sr2MgOs6_noSOC.struct`. We use the following:

.. literalinclude:: images_scripts/Sr2MgOsO6_noSOC.struct

The DFT calculation is done as usual, for instance you can use for the initialisation


   init -b -vxc 5 -numk 2000 

This is setting up a non-magnetic calculation, using the LDA and 2000 k-points in the full Brillouin zone. As usual, we start the DFT self consistent cycle by the Wien2k script ::

  run

Wannier orbitals
----------------

As a next step, we calculate localised orbitals for the t2g orbitals. We use the same input file for :program:`dmftproj` as it was used in the :ref:`documentation`:

.. literalinclude:: images_scripts/Sr2MgOsO6_noSOC.indmftpr

Note that, due to the distortions in the crystal structure, we need to include all five d orbitals in the calculation (line 8 in the input file above).

To prepare the input data for :program:`dmftproj` we execute lapw2 with the `-almd` option ::
   
   x lapw2 -almd 

Then  :program:`dmftproj` is executed in its default mode (i.e. without spin-polarization or spin-orbit included) ::

   dmftproj 

This program produces the necessary files for the conversion to the hdf5 file structure. This is done using
the python module :class:`Wien2kConverter <dft.converters.wien2k_converter.Wien2kConverter>`. A simple python script that initialises the converter is::

  from triqs_dft_tools.converters.wien2k_converter import *
  Converter = Wien2kConverter(filename = "Sr2MgOsO6_noSOC")

After initializing the interface module, we can now convert the input
text files to the hdf5 archive by::

  Converter.convert_dft_input()

This reads all the data, and stores everything that is necessary for the DMFT calculation in the file :file:`Sr2MgOsO6_noSOC.h5`.

[CONTINUE HERE]

The DMFT calculation
====================

