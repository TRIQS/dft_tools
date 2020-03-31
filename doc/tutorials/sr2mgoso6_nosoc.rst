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

After the SC cycled finished, you can calculate the DOS. It should look like what you can see in this figure:

.. image:: images_scripts/Sr2MgOsO6_noSOC_DOS.png
    :width: 700
    :align: center


Wannier orbitals
----------------

As a next step, we calculate localised orbitals for the d orbitals. We use this input file for :program:`dmftproj`:

.. literalinclude:: images_scripts/Sr2MgOsO6_noSOC.indmftpr

Note that, due to the distortions in the crystal structure, we need to include all five d orbitals in the calculation (line 8 in the input file above).

To prepare the input data for :program:`dmftproj` we execute lapw2 with the `-almd` option ::
   
   x lapw2 -almd 

Then  :program:`dmftproj` is executed in its default mode (i.e. without spin-polarization or spin-orbit included) ::

   dmftproj 

At the end of the run you see the density matrix in Wannier space:

.. literalinclude:: images_scripts/Sr2MgOsO6_noSOC.outdmftpr

As you can see, there are off-diagonal elements between the :math:`d_{x^2-y^2}` and the :math:`d_{xy}` orbital.

We convert the output to the hdf5 archive, using 
the python module :class:`Wien2kConverter <dft.converters.wien2k_converter.Wien2kConverter>`. A simple python script doing this is::

  from triqs_dft_tools.converters.wien2k_converter import *
  Converter = Wien2kConverter(filename = "Sr2MgOsO6_noSOC")
  Converter.convert_dft_input()

This reads all the data, and stores everything that is necessary for the DMFT calculation in the file :file:`Sr2MgOsO6_noSOC.h5`.


The DMFT calculation
====================

Before starting the DMFT calculation it is beneficial to look a bit more closely at the block structure of the problem. Eventually, we want to use a basis that is as diagonal as possible, and include only the partially filled orbitals in the correlated problem. All this can be done using the functionalities of the :class:`BlockStructure <dft.block_structure.BlockStructure>` class, see section :ref:`blockstructure`.

We first initialise the SumK class::

    from triqs_dft_tools.sumk_dft import *
    SK = SumkDFT(hdf_file='Sr2MgOsO6_noSOC.h5',use_dft_blocks=True)

The flag *use_dft_blocks=True* determines, as usual, the smallest possible blocks with non-zero entries, and initialises them as *solver* block structure. In order to disentangle the :math:`d_{x^2-y^2}` and the :math:`d_{xy}` orbitals, and pick out only the partially filled one, we do a transformation into a basis where the local Hamiltonian is diagonal::

    mat = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal',calc_in_solver_blocks=True)

We can look at the diagonalisation matrix, it is::

    >>> print mat[0]['down']
    [[ 1.   +0.j     0.   +0.j     0.   +0.j     0.   +0.j     0.   +0.j   ]
     [ 0.   +0.j    -0.985+0.j     0.   -0.173j  0.   +0.j     0.   +0.j   ]
     [ 0.   +0.j     0.173+0.j     0.   -0.985j  0.   +0.j     0.   +0.j   ]
     [ 0.   +0.j     0.   +0.j     0.   +0.j     1.   +0.j     0.   +0.j   ]
     [ 0.   +0.j     0.   +0.j     0.   +0.j     0.   +0.j     1.   +0.j   ]]
    >>> 

This transformation is already stored in the SK.block_structure class. The next step is actually not needed for a DMFT calculation, but lets see what the transformation does to the local Hamiltonian. We can calculate it before rotation, rotate it, and look at the 2x2 block with off-diagonals::

    eal  = SK.eff_atomnic_levels()
    eal2 = SK.block_structure.convert_matrix(eal[0],space_from='sumk', space_to='solver')

    print eal[0]['up'][1:3,1:3]                # prints the 2x2 block with offiagonals
    [[ 0.391+0.j    -0.   -0.815j]
     [-0.   +0.815j  4.884-0.j   ]]

    print eal2['up_1']                         # prints the 2x2 block after rotation
    [[0.247-0.j 0.   +0.j]
     [0.   -0.j 5.028+0.j]]

So the local Hamiltonian has been diagonalised. From the other entries we can see that the *up_0* block and the [1,1] entry of the *up_1* block correspond to :math:`e_g`-like orbitals, and the others are the 
:math:`t_{2g}` orbitals that we want to keep. 

[TO BE CONTINUED...]
