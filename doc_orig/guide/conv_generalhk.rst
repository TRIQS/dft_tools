.. _convgeneralhk:

A general H(k)
==============

In addition to the more extensive Wien2k, VASP, and W90 converters,
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

  from triqs_dft_tools.converters.hk import *
  Converter = HkConverter(filename = hkinputfile)
  Converter.convert_dft_input()

where :file:`hkinputfile` is the name of the input file described
above. This produces the hdf file that you need for a DMFT calculation.

For more options of this converter, have a look at the
:ref:`refconverters` section of the reference manual.


