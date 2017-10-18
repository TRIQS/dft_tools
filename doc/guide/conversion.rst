.. _conversion:

Orbital construction and conversion
===================================

The first step for a DMFT calculation is to provide the necessary
input based on a DFT calculation. We will not review how to do the DFT
calculation here in this documentation, but refer the user to the
documentation and tutorials that come with the actual DFT
package. Here, we will describe how to use output created by Wien2k,
as well as how to use the light-weight general interface.

Interface with Wien2k
---------------------

We assume that the user has obtained a self-consistent solution of the
Kohn-Sham equations. We further have to require that the user is
familiar with the main in/output files of Wien2k, and how to run
the DFT code.

Conversion for the DMFT self-consistency cycle
""""""""""""""""""""""""""""""""""""""""""""""

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

Let us take the compound SrVO3, a commonly used
example for DFT+DMFT calculations. The input file for
:program:`dmftproj` looks like 

.. literalinclude:: images_scripts/SrVO3.indmftpr

The first three lines give the number of inequivalent sites, their
multiplicity (to be in accordance with the Wien2k *struct* file) and
the maximum orbital quantum number :math:`l_{max}`. In our case our
struct file contains the atoms in the order Sr, V, O.

Next we have to
specify for each of the inequivalent sites, whether we want to treat
their orbitals as correlated or not. This information is given by the
following 3 to 5 lines:

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

The last line gives the energy window, relative to the Fermi energy,
that is used for the projective Wannier functions. Note that, in
accordance with Wien2k, we give energies in Rydberg units!

After setting up this input file, you run:

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
use the python module :class:`Wien2kConverter <dft.converters.wien2k_converter.Wien2kConverter>`. It is initialized as::

  from pytriqs.applications.dft.converters.wien2k_converter import *
  Converter = Wien2kConverter(filename = case)

The only necessary parameter to this construction is the parameter `filename`.
It has to be the root of the files produces by dmftproj. For our
example, the :program:`Wien2k` naming convention is that all files are
called the same, for instance
:file:`SrVO3.*`, so you would give `filename = "SrVO3"`. The constructor opens
an hdf5 archive, named :file:`case.h5`, where all the data is
stored. For other parameters of the constructor please visit the
:ref:`refconverters` section of the reference manual.

After initializing the interface module, we can now convert the input
text files to the hdf5 archive by::

  Converter.convert_dft_input()

This reads all the data, and stores it in the file :file:`case.h5`. 
In this step, the files :file:`case.ctqmcout` and
:file:`case.symqmc` 
have to be present in the working directory.

After this step, all the necessary information for the DMFT loop is
stored in the hdf5 archive, where the string variable
`Converter.hdf_filename` gives the file name of the archive.

At this point you should use the method :meth:`dos_wannier_basis <dft.sumk_dft_tools.SumkDFTTools.dos_wannier_basis>`
contained in the module :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>` to check the density of
states of the Wannier orbitals (see :ref:`analysis`).

You have now everything for performing a DMFT calculation, and you can
proceed with the section on :ref:`single-shot DFT+DMFT calculations <singleshot>`.

Data for post-processing
""""""""""""""""""""""""

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
"""""""""""""""""""""""""""""""

For the transport calculations, the situation is a bit more involved,
since we need also the :program:`optics` package of Wien2k. Please
look at the section on :ref:`Transport` to see how to do the necessary
steps, including the conversion.
  
Interface with VASP
---------------------

.. warning::
  The VASP interface is in the alpha-version and the VASP part of it is not
  yet publicly released. The documentation may, thus, be subject to changes
  before the final release.

The interface with VASP relies on new options introduced since
version 5.4.x. The output of raw (non-normalized) projectors is
controlled by an INCAR option LOCPROJ whose complete syntax is described in
VASP documentaion.

The definition of a projector set starts with specifying which sites
and which local states we are going to project onto.
This information is provided by option LOCPROJ

  | `LOCPROJ = <sites> : <shells> : <projector type>`

where `<sites>` represents a list of site indices separated by spaces,
with the indices corresponding to the site position in the POSCAR file;
`<shells>` specifies local states (e.g. :math:`s`, :math:`p`, :math:`d`,
:math:`d_{x^2-y^2}`, etc.);
`<projector type>` chooses a particular type of the local basis function.

Some projector types also require parameters `EMIN`, `EMAX` in INCAR to
be set to define an (approximate) energy window corresponding to the 
valence states.

When either a self-consistent (`ICHARG < 10`) or a non-self-consistent
(`ICHARG >= 10`) calculation is done VASP produces file `LOCPROJ` which
will serve as the main input for the conversion routine.


Conversion for the DMFT self-consistency cycle
""""""""""""""""""""""""""""""""""""""""""""""

In order to use the projectors generated by VASP for defining an
impurity problem they must be processed, i.e. normalized, possibly
transformed, and then converted to a format suitable for DFT_tools scripts.

The processing of projectors is performed by the program :program:`plovasp`
invoked as

  | `plovasp <plo.cfg>`

where `<plo.cfg>` is a input file controlling the conversion of projectors.

The format of input file `<plo.cfg>` is described in details in
:ref:`plovasp`. Here we just give a simple example for the case
of SrVO3:

.. literalinclude:: images_scripts/srvo3.cfg

A projector shell is defined by a section `[Shell 1]` where the number
can be arbitrary and used only for user convenience. Several
parameters are required

- **IONS**: list of site indices which must be a subset of indices
  given earlier in `LOCPROJ`.
- **LSHELL**: :math:`l`-quantum number of the projector shell; the corresponding
  orbitals must be present in `LOCPROJ`.
- **EWINDOW**: energy window in which the projectors are normalized;
  note that the energies are defined with respect to the Fermi level.
  
Option **TRANSFORM** is optional but here it is specified to extract
only three :math:`t_{2g}` orbitals out of five `d` orbitals given by 
:math:`l = 2`.

 
A general H(k)
--------------

In addition to the more complicated Wien2k converter,
:program:`DFTTools` contains also a light converter. It takes only
one inputfile, and creates the necessary hdf outputfile for
the DMFT calculation. The header of this input file has a defined
format, an example is the following (do not use the text/comments in your
input file):

.. literalinclude:: images_scripts/case.hk

The lines of this header define
		    
#. Number of :math:`\mathbf{k}`-points used in the calculation
#. Electron density for setting the chemical potential
#. Number of total atomic shells in the hamiltonian matrix. In short,
   this gives the number of lines described in the following. IN the
   example file give above this number is 2.
#. The next line(s) contain four numbers each: index of the atom, index
   of the equivalent shell, :math:`l` quantum number, dimension
   of this shell. Repeat this line for each atomic shell, the number
   of the shells is given in the previous line.

   In the example input file given above, we have two inequivalent
   atomic shells, one on atom number 1 with a full d-shell (dimension 5),
   and one on atom number 2 with one p-shell (dimension 3).

   Other examples for these lines are:

   #. Full d-shell in a material with only one correlated atom in the
      unit cell (e.g. SrVO3). One line is sufficient and the numbers
      are `1 1 2 5`.
   #. Full d-shell in a material with two equivalent atoms in the unit
      cell (e.g. FeSe): You need two lines, one for each equivalent
      atom. First line is `1 1 2 5`, and the second line is
      `2 1 2 5`. The only difference is the first number, which tells on
      which atom the shell is located. The second number is the
      same in both lines, meaning that both atoms are equivalent.
   #. t2g orbitals on two non-equivalent atoms in the unit cell: Two
      lines again. First line is `1 1 2 3`, second line `2 2 2 3`. The
      difference to the case above is that now also the second number
      differs. Therefore, the two shells are treated independently in
      the calculation.
   #. d-p Hamiltonian in a system with two equivalent atoms each in
      the unit cell (e.g. FeSe has two Fe and two Se in the unit
      cell). You need for lines. First line `1 1 2 5`, second
      line
      `2 1 2 5`. These two lines specify Fe as in the case above. For the p
      orbitals you need line three as `3 2 1 3` and line four
      as `4 2 1 3`. We have 4 atoms, since the first number runs from 1 to 4,
      but only two inequivalent atoms, since the second number runs
      only form 1 to 2.
   
   Note that the total dimension of the hamiltonian matrices that are
   read in is the sum of all shell dimensions that you specified. For
   example number 4 given above we have a dimension of 5+5+3+3=16. It is important
   that the order of the shells that you give here must be the same as
   the order of the orbitals in the hamiltonian matrix. In the last
   example case above the code assumes that matrix index 1 to 5
   belongs to the first d shell, 6 to 10 to the second, 11 to 13 to
   the first p shell, and 14 to 16 the second p shell. 
   
#. Number of correlated shells in the hamiltonian matrix, in the same
   spirit as line 3.

#. The next line(s) contain six numbers: index of the atom, index
   of the equivalent shell, :math:`l` quantum number, dimension
   of the correlated shells, a spin-orbit parameter, and another
   parameter defining interactions. Note that the latter two
   parameters are not used at the moment in the code, and only kept
   for compatibility reasons. In our example file we use only the
   d-shell as correlated, that is why we have only one line here.
      
#. The last line contains several numbers: the number of irreducible
   representations, and then the dimensions of the irreps. One
   possibility is as the example above, another one would be 2
   2 3. This would mean, 2 irreps (eg and t2g), of dimension 2 and 3,
   resp.

After these header lines, the file has to contain the Hamiltonian
matrix in orbital space. The standard convention is that you give for
each :math:`\mathbf{k}`-point first the matrix of the real part, then the
matrix of the imaginary part, and then move on to the next :math:`\mathbf{k}`-point.

The converter itself is used as::

  from pytriqs.applications.dft.converters.hk_converter import *
  Converter = HkConverter(filename = hkinputfile)
  Converter.convert_dft_input()
  
where :file:`hkinputfile` is the name of the input file described
above. This produces the hdf file that you need for a DMFT calculation.

For more options of this converter, have a look at the
:ref:`refconverters` section of the reference manual.


Wannier90 Converter
-------------------

Using this converter it is possible to convert the output of
`wannier90 <http://wannier.org>`_
Maximally Localized Wannier Functions (MLWF) and create a HDF5 archive
suitable for one-shot DMFT calculations with the
:class:`SumkDFT <dft.sumk_dft.SumkDFT>` class.

The user must supply two files in order to run the Wannier90 Converter:

#. The file :file:`seedname_hr.dat`, which contains the DFT Hamiltonian
   in the MLWF basis calculated through :program:`wannier90` with ``hr_plot = true``
   (please refer to the :program:`wannier90` documentation).
#. A file named :file:`seedname.inp`, which contains the required
   information about the :math:`\mathbf{k}`-point mesh, the electron density,
   the correlated shell structure, ... (see below).

Here and in the following, the keyword ``seedname`` should always be intended
as a placeholder for the actual prefix chosen by the user when creating the
input for :program:`wannier90`.
Once these two files are available, one can use the converter as follows::

    from pytriqs.applications.dft.converters import Wannier90Converter
    Converter = Wannier90Converter(seedname='seedname')
    Converter.convert_dft_input()

The converter input :file:`seedname.inp` is a simple text file with
the following format (do not use the text/comments in your input file):

.. literalinclude:: images_scripts/LaVO3_w90.inp

The example shows the input for the perovskite crystal of LaVO\ :sub:`3`
in the room-temperature `Pnma` symmetry. The unit cell contains four
symmetry-equivalent correlated sites (the V atoms) and the total number
of electrons per unit cell is 8 (see second line).
The first line specifies how to generate the :math:`\mathbf{k}`-point
mesh that will be used to obtain :math:`H(\mathbf{k})`
by Fourier transforming :math:`H(\mathbf{R})`.
Currently implemented options are:

* :math:`\Gamma`-centered uniform grid with dimensions
  :math:`n_{k_x} \times n_{k_y} \times n_{k_z}`;
  specify ``0`` followed by the three grid dimensions,
  like in the example above
* :math:`\Gamma`-centered uniform grid with dimensions 
  automatically determined by the converter (from the number of
  :math:`\mathbf{R}` vectors found in :file:`seedname_hr.dat`);
  just specify ``-1``

Inside :file:`seedname.inp`, it is crucial to correctly specify the
correlated shell structure, which depends on the contents of the
:program:`wannier90` output :file:`seedname_hr.dat` and on the order
of the MLWFs contained in it. In this example we have four lines for the
four V atoms. The MLWFs were constructed for the t\ :sub:`2g` subspace, and thus
we set ``l`` to 2 and ``dim`` to 3 for all V atoms. Further the spin-orbit coupling (``SO``)
is set to 0 and ``irep`` to 0.
As in this example all 4 V atoms are equivalent we set ``sort`` to 0. We note
that, e.g., for a magnetic DMFT calculation the correlated atoms can be made
inequivalent at this point by using different values for ``sort``.

The number of MLWFs must be equal to, or greater than the total number
of correlated orbitals (i.e., the sum of all ``dim`` in :file:`seedname.inp`).
If the converter finds fewer MLWFs inside :file:`seedname_hr.dat`, then it
stops with an error; if it finds more MLWFs, then it assumes that the
additional MLWFs correspond to uncorrelated orbitals (e.g., the O-\ `2p` shells).
When reading the hoppings :math:`\langle w_i | H(\mathbf{R}) | w_j \rangle`
(where :math:`w_i` is the :math:`i`-th MLWF), the converter also assumes that
the first indices correspond to the correlated shells (in our example,
the V-t\ :sub:`2g` shells). Therefore, the MLWFs corresponding to the
uncorrelated shells (if present) must be listed **after** those of the 
correlated shells.
With the :program:`wannier90` code, this can be achieved by listing the
projections for the uncorrelated shells after those for the correlated shells.
In our `Pnma`-LaVO\ :sub:`3` example, for instance, we could use::

    Begin Projections
     V:l=2,mr=2,3,5:z=0,0,1:x=-1,1,0
     O:l=1:mr=1,2,3:z=0,0,1:x=-1,1,0
    End Projections

where the ``x=-1,1,0`` option indicates that the V--O bonds in the octahedra are
rotated by (approximatively) 45 degrees with respect to the axes of the `Pbnm` cell.

The converter will analyze the matrix elements of the local Hamiltonian
to find the symmetry matrices `rot_mat` needed for the global-to-local
transformation of the basis set for correlated orbitals
(see section :ref:`hdfstructure`).
The matrices are obtained by finding the unitary transformations that diagonalize
:math:`\langle w_i | H_I(\mathbf{R}=0,0,0) | w_j \rangle`, where :math:`I` runs
over the correlated shells and `i,j` belong to the same shell (more details elsewhere...).
If two correlated shells are defined as equivalent in :file:`seedname.inp`,
then the corresponding eigenvalues have to match within a threshold of 10\ :sup:`-5`,
otherwise the converter will produce an error/warning.
If this happens, please carefully check your data in :file:`seedname_hr.dat`.
This method might fail in non-trivial cases (i.e., more than one correlated
shell is present) when there are some degenerate eigenvalues:
so far tests have not shown any issue, but one must be careful in those cases
(the converter will print a warning message).

The current implementation of the Wannier90 Converter has some limitations:

* Since :program:`wannier90` does not make use of symmetries (symmetry-reduction
  of the :math:`\mathbf{k}`-point grid is not possible), the converter always
  sets ``symm_op=0`` (see the :ref:`hdfstructure` section).
* No charge self-consistency possible at the moment.
* Calculations with spin-orbit (``SO=1``) are not supported.
* The spin-polarized case (``SP=1``) is not yet tested.
* The post-processing routines in the module
  :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`
  were not tested with this converter.
* ``proj_mat_all`` are not used, so there are no projectors onto the
  uncorrelated orbitals for now.

     
MPI issues
----------

The interface packages are written such that all the file operations
are done only on the master node. In general, the philosophy of the
package is that whenever you read in something from the archive
yourself, you have to *manually* broadcast it to the nodes. An
exception to this rule is when you use routines from :class:`SumkDFT <dft.sumk_dft.SumkDFT>`
or :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`, where the broadcasting is done for you.

Interfaces to other packages
----------------------------

Because of the modular structure, it is straight forward to extend the :ref:`TRIQS <triqslibs:welcome>` package
in order to work with other band-structure codes. The only necessary requirement is that
the interface module produces an hdf5 archive, that stores all the data in the specified
form. For the details of what data is stored in detail, see the
:ref:`hdfstructure` part of the reference manual.
