
The interface
=============


The basic function of the interface to the Wien2k program package is
to take the output of the program that constructs the projected local
orbitals (:program:`dmftproj`, for documentation see :download:`TutorialDmftproj.pdf <TutorialDmftproj.pdf>`), and to store all the necessary information into
an hdf5 file. This latter file is then used to do the DMFT calculation. The
reason for this structure is that this enables the user to have everything
that is necessary to reproduce the calculation in one single hdf5 arxive.

.. index:: Interface to Wien2k

.. _interfacetowien:

The interface to Wien2k
-----------------------

As explained above, this interface produces an hdf5 arxive out of the files that
were written by the band structure package :program:`Wien2k/dmftproj`. 
For this purpose we
use the python module :class:`Wien2kConverter`. It is initialised as::

  from pytriqs.applications.dft.converters.wien2k_converter import *
  Converter = Wien2kConverter(filename = material_of_interest)

The only necessary parameter to this construction is the parameter `filename`.
It has to be the root of the files produces by dmftproj. For example, if you did a 
calculation for TiO, the :program:`Wien2k` naming convention is that all files are called 
:file:`TiO.*`, so you would give `filename = "TiO"`. The constructor opens
an hdf5 arxive, named :file:`material_of_interest.h5`, where all the data is stored.

There are three optional parameters to the Constructor:

  * `lda_subgrp`: We store all data in sub groups of the hdf5 arxive. For the main data
    that is needed for the DMFT loop, we use the sub group specified by this optional parameter.
    If it is not given, the standard value `SumK_LDA` is used as sub group name.
  * `symm_subgrp`: In this sub group we store all the data for applying the symmetry 
    operations in the DMFT loop. Standard value is `SymmCorr`.
  * `repacking`: If true, and the hdf5 file already exists, the system command :program:`h5repack` 
    is invoked. This command ensures a minimal file size of the hdf5
    file. Standard value is `False`. If you want to use this, be sure
    that :program:`h5repack` is in your path variable!

After initialising the interface module, we can now convert the input text files into the
hdf5 arxive by::

  Converter.convert_dmft_input()

This reads all the data, and stores it in the sub group `lda_subgrp`, as discussed above. 
In this step, the files :file:`material_of_interest.ctqmcout` and :file:`material_of_interest.symqmc`
have to be present in the working directory.

After this step, all the necessary information for the DMFT loop is stored in the hdf5 arxive, where
the string variable `Converter.hdf_file` gives the file name of the arxive.
You can now proceed with :ref:`LDADMFTmain`.


Data for post-processing
------------------------

In order to calculate some properties using the DMFT self energy, several other routines are
used in order to convert the necessary input from :program:`Wien2k/dmftproj`. For instance, for 
calculating the partial density of states or partial charges consistent with the definition
of :program:`Wien2k`, you have to use::

  Converter.convert_parproj_input()

This reads the files :file:`material_of_interest.parproj` and :file:`material_of_interest.sympar`.
Again, there are two optional parameters

  * `par_proj_subgrp`: The sub group, where the data for the partial projectors is stored. Standard
    is `SumK_LDA_ParProj`.
  * `symm_par_subgrp`: Sub group for the symmetry operations, standard value is `SymmPar`.

Another routine of the class allows to read the input for plotting the momentum-resolved
spectral function. It is done by::
  
  Converter.convert_bands_input()

The optional parameter, which tells the routine where to store the data is here `bands_subgrp`, 
and its standard value is `SumK_LDA_Bands`.

After having converted this input, you can further proceed with the :ref:`analysis`.

MPI issues
----------

The interface package is written such that all the operations are done only on the master node.
The broadcasting to the nodes has to be done by hand. The :class:`SumkLDA`, described in the
following section, takes care of this automatically.

Interfaces to other packages
----------------------------

Because of the modular structure, it is straight forward to extend the TRIQS package 
in order to work with other band-structure codes. The only necessary requirement is that 
the interface module produces an hdf5 arxive, that stores all the data in the specified
form. For the details of what data is stored in detail, see the reference manual.
