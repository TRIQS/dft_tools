.. _convWien2k:

Interface with Wien2k
=====================

We assume that the user has obtained a self-consistent solution of the
Kohn-Sham equations. We further have to require that the user is
familiar with the main in/output files of Wien2k, and how to run
the DFT code.

Conversion for the DMFT self-consistency cycle
----------------------------------------------

First, we have to write the necessary
quantities into a file that can be processed further by invoking in a
shell the command

  `x lapw2 -almd`

We note that any other flag for lapw2, such as -c or -so (for
spin-orbit coupling) has to be added also to this line. This creates
some files that we need for the Wannier orbital construction.

The orbital construction itself is done by the Fortran program
:program:`dmftproj`. For an extensive manual to this program see
:download:`TutorialDmftproj.pdf <images_scripts/TutorialDmftproj.pdf>`.
Here we will only describe the basic steps.

In the following, we use SrVO3 as an example to explain the
input file :file:`case.indmftpr` for :program:`dmftproj`.
A full tutorial on SrVO3 is available in the :ref:`SrVO3 tutorial <SrVO3>`.

.. literalinclude:: ../tutorials/images_scripts/SrVO3.indmftpr

The first three lines give the number of inequivalent sites, their
multiplicity (to be in accordance with the Wien2k *struct* file) and
the maximum orbital quantum number :math:`l_{max}`. In our case our
struct file contains the atoms in the order Sr, V, O.

Next we have to specify for each of the inequivalent sites, whether
we want to treat their orbitals as correlated or not. This information
is given by the following 3 to 5 lines:

#. We specify which basis set is used (complex or cubic
   harmonics).
#. The four numbers refer to *s*, *p*, *d*, and *f* electrons,
   resp. Putting 0 means doing nothing, putting 1 will calculate
   **unnormalized** projectors in compliance with the Wien2k
   definition. The important flag is 2, this means to include these
   electrons as correlated electrons, and calculate normalized Wannier
   functions for them. In the example above, you see that only for the
   vanadium *d* we set the flag to 2. If you want to do simply a DMFT
   calculation, then set everything to 0, except one flag 2 for the
   correlated electrons.
#. In case you have a irrep splitting of the correlated shell, you can
   specify here how many irreps you have. You see that we put 2, since
   eg and t2g symmetries are irreps in this cubic case. If you don't
   want to use this splitting, just put 0.
#. (optional) If you specifies a number different from 0 in above line, you have
   to tell now, which of the irreps you want to be treated
   correlated. We want to t2g, and not the eg, so we set 0 for eg and
   1 for t2g. Note that the example above is what you need in 99% of
   the cases when you want to treat only t2g electrons. For eg's only
   (e.g. nickelates), you set 10 and 01 in this line.
#. (optional) If you have specified a correlated shell for this atom,
   you have to tell if spin-orbit coupling should be taken into
   account. 0 means no, 1 is yes.

These lines have to be repeated for each inequivalent atom.

The last line gives the lower and upper limit of the energy window,
relative to the Fermi energy, which is used for the projective Wannier functions.
Note that, in accordance with Wien2k, we give energies in Rydberg units!

The third number is an optional flag to switch between different modes:

#. 0: The projectors are constructed for the given energy window. The number
   of bands within the window is usually different at each k-point which
   will be reflected by the projectors, too. This is the default mode
   which is also used if no mode flag is provided.
#. 1: The lowest and highest band indices within the given energy window
   are calculated. The resulting indices are used at all k-points.
   Bands which fall within the window only in some parts of the Brillouin zone
   are fully taken into account. Keep in mind that a different set of k-points
   or the -band option can change the lower or upper index. This can be avoided
   by using mode 2.
#. 2: In this mode the first two values of the line are interpreted as lower
   and upper band indices to be included in the projective subspace. For example,
   if the line reads `21 23 2`, bands number 21, 22 and 23 are included at all
   k-points. For SrVO3 this corresponds to the t2g bands around the Fermi energy.
   The lowest possible index is 1. Note that the indices need to be provided as integer.

In all modes the used energy range, i.e. band range, is printed to the
:program:`dmftproj` output.

After setting up the :file:`case.indmftpr` input file, you run:

  `dmftproj`

Again, adding possible flags like -so for spin-orbit coupling. This
program produces the following files (in the following, take *case* as
the standard Wien2k place holder, to be replaced by the actual working
directory name):

 * :file:`case.ctqmcout` and :file:`case.symqmc` containing projector
   operators and symmetry operations for orthonormalized Wannier
   orbitals, respectively.
 * :file:`case.parproj` and :file:`case.sympar` containing projector
   operators and symmetry operations for uncorrelated states,
   respectively. These files are needed for projected
   density-of-states or spectral-function calculations in
   post-processing only.
 * :file:`case.oubwin` needed for the charge density recalculation in
   the case of fully self-consistent DFT+DMFT run (see below).

Now we convert these files into an hdf5 file that can be used for the
DMFT calculations. For this purpose we
use the python module :class:`Wien2kConverter <dft.converters.wien2k.Wien2kConverter>`. It is initialized as::

  from triqs_dft_tools.converters.wien2k import *
  Converter = Wien2kConverter(filename = case)

The only necessary parameter to this construction is the parameter `filename`.
It has to be the root of the files produces by dmftproj. For our
example, the :program:`Wien2k` naming convention is that all files have the
same name, but different extensions, :file:`case.*`. The constructor opens
an hdf5 archive, named :file:`case.h5`, where all relevant data will be
stored. For other parameters of the constructor please visit the
:ref:`refconverters` section of the reference manual.

After initializing the interface module, we can now convert the input
text files to the hdf5 archive by::

  Converter.convert_dft_input()

This reads all the data, and stores it in the file :file:`case.h5`.
In this step, the files :file:`case.ctqmcout` and
:file:`case.symqmc` have to be present in the working directory.

After this step, all the necessary information for the DMFT loop is
stored in the hdf5 archive, where the string variable
`Converter.hdf_filename` gives the file name of the archive.

At this point you should use the method :meth:`dos_wannier_basis <dft.sumk_dft_tools.SumkDFTTools.dos_wannier_basis>`
contained in the module :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>` to check the density of
states of the Wannier orbitals (see :ref:`analysis`).

You have now everything for performing a DMFT calculation, and you can
proceed with the section on :ref:`single-shot DFT+DMFT calculations <singleshot>`.

Data for post-processing
------------------------

In case you want to do post-processing of your data using the module
:class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`, some more files
have to be converted to the hdf5 archive. For instance, for
calculating the partial density of states or partial charges
consistent with the definition of :program:`Wien2k`, you have to invoke::

  Converter.convert_parproj_input()

This reads and converts the files :file:`case.parproj` and
:file:`case.sympar`.

If you want to plot band structures, one has to do the
following. First, one has to do the Wien2k calculation on the given
:math:`\mathbf{k}`-path, and run :program:`dmftproj` on that path:

  |  `x lapw1 -band`
  |  `x lapw2 -band -almd`
  |  `dmftproj -band`


Again, maybe with the optional additional extra flags according to
Wien2k. Now we use a routine of the converter module allows to read
and convert the input for :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`::

  Converter.convert_bands_input()

After having converted this input, you can further proceed with the
:ref:`analysis`. For more options on the converter module, please have
a look at the :ref:`refconverters` section of the reference manual.

Data for transport calculations
-------------------------------

For the transport calculations, the situation is a bit more involved,
since we need also the :program:`optics` package of Wien2k. Please
look at the section on :ref:`Transport` to see how to do the necessary
steps, including the conversion.


