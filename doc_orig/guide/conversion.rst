.. _conversion:

Supported interfaces
====================

The first step for a DMFT calculation is to provide the necessary
input based on a DFT calculation. We will not review how to do the DFT
calculation here in this documentation, but refer the user to the
documentation and tutorials that come with the actual DFT
package. At the moment, there are two full charge self consistent interfaces, for the
Wien2k and the VASP DFT packages, resp. In addition, there is an interface to Wannier90, as well
as a light-weight general-purpose interface. In the following, we will describe the usage of these
conversion tools.

.. toctree::
   :maxdepth: 3

   conv_wien2k
   conv_vasp
   conv_W90
   conv_generalhk

MPI issues
==========

The interface packages are written such that all the file operations
are done only on the master node. In general, the philosophy of the
package is that whenever you read in something from the archive
yourself, you have to *manually* broadcast it to the nodes. An
exception to this rule is when you use routines from :class:`SumkDFT <dft.sumk_dft.SumkDFT>`
or :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`, where the broadcasting is done for you.


Interfaces to other packages
============================

Because of the modular structure, it is straight forward to extend the :ref:`TRIQS <triqslibs:welcome>` package
in order to work with other band-structure codes. The only necessary requirement is that
the interface module produces an hdf5 archive, that stores all the data in the specified
form. For the details of what data is stored in detail, see the
:ref:`hdfstructure` part of the reference manual.
