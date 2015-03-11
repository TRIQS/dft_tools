
The interface
=============

The basic function of the interface to the Wien2k program package is to take
the output of the program that constructs the projected local orbitals
(:program:`dmftproj`, for documentation see 
:download:`TutorialDmftproj.pdf <TutorialDmftproj.pdf>`), 
and to store all the necessary information into an hdf5 file. This latter file
is then used to do the DMFT calculation. The reason for this structure is that
this enables the user to have everything that is necessary to reproduce the
calculation in one single hdf5 archive.

.. index:: Interface to Wien2k

.. _interfacetowien:

The interface to Wien2k
-----------------------

As explained above, this interface produces an hdf5 archive out of the files that
were written by the band structure package :program:`Wien2k/dmftproj`. 
For this purpose we
use the python module :class:`Wien2kConverter`. It is initialised as::

  from pytriqs.applications.dft.converters.wien2k_converter import *
  Converter = Wien2kConverter(filename = material_of_interest)

The only necessary parameter to this construction is the parameter `filename`.
It has to be the root of the files produces by dmftproj. For example, if you did a 
calculation for TiO, the :program:`Wien2k` naming convention is that all files are called 
:file:`TiO.*`, so you would give `filename = "TiO"`. The constructor opens
an hdf5 archive, named :file:`material_of_interest.h5`, where all the data is stored.

These are the parameters to the Constructor:

=========================   ============================  ===========================================================================
Name                        Type, Default                 Meaning
=========================   ============================  ===========================================================================
filename                    String                        Material being studied, corresponding to the :program:`Wien2k` file names.
                                                          The constructor stores the data in the hdf5 archive :file:`material_of_interest.h5`.
dft_subgrp                  String, dft_input             hdf5 subgroup containing required DFT data
symmcorr_subgrp             String, dft_symmcorr_input    hdf5 subgroup containing all necessary data to apply
                                                          the symmetry operations in the DMFT loop
repacking                   Boolean, False                Does the hdf5 file already exist and should the :program:`h5repack` be 
                                                          invoked to ensures a minimal archive file size? 
                                                          Note that the :program:`h5repack` must be in your path variable!
=========================   ============================  ===========================================================================

After initialising the interface module, we can now convert the input text files into the
hdf5 archive by::

  Converter.convert_dft_input()

This reads all the data, and stores it in the subgroup `dft_subgrp`, as discussed above. 
In this step, the files :file:`material_of_interest.ctqmcout` and :file:`material_of_interest.symqmc`
have to be present in the working directory.

After this step, all the necessary information for the DMFT loop is stored in the hdf5 archive, where
the string variable `Converter.hdf_file` gives the file name of the archive.
You can now proceed with :ref:`DFTDMFTmain`.


Data for post-processing
------------------------

In order to calculate some properties using the DMFT self energy, several other routines are
used in order to convert the necessary input from :program:`Wien2k/dmftproj`. For instance, for 
calculating the partial density of states or partial charges consistent with the definition
of :program:`Wien2k`, you have to use::

  Converter.convert_parproj_input()

This reads the files :file:`material_of_interest.parproj` and :file:`material_of_interest.sympar`.
Again, there are two optional parameters

=========================   ============================  ===========================================================================
Name                        Type, Default                 Meaning
=========================   ============================  ===========================================================================
parproj_subgrp              String, dft_parproj_input     hdf5 subgroup containing partial projectors data.
symmpar_subgrp              String, dft_symmpar_input     hdf5 subgroup containing symmetry operations data.
=========================   ============================  ===========================================================================

Another routine of the class allows to read the input for plotting the momentum-resolved
spectral function. It is done by::
  
  Converter.convert_bands_input()

The optional parameter that controls where the data is stored is `bands_subgrp`, 
with the default value `dft_bands_input`. Note however that you need to run "dmftproj -band" to produce the
necessary outband file. The casename.indmftpr file needs an additional line with E_fermi 
(obtainable from casename.qtl).

After having converted this input, you can further proceed with the :ref:`analysis`.

MPI issues
----------

The interface package is written such that all the operations are done only on the master node.
The broadcasting to the nodes has to be done by hand. The :class:`SumkDFT`, described in the
following section, takes care of this automatically.

Interfaces to other packages
----------------------------

Because of the modular structure, it is straight forward to extend the TRIQS package 
in order to work with other band-structure codes. The only necessary requirement is that 
the interface module produces an hdf5 archive, that stores all the data in the specified
form. For the details of what data is stored in detail, see the reference manual.
