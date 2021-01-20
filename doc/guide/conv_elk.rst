.. _convElk:

Interface with Elk
=====================

This is the first iteration of the Elk-TRIQS interface, so certain inputs may change in later updates. The Elk part of the interface is not currently in the main distribution, but it can be found :ref:`here <https://github.com/AlynJ/Elk_interface-TRIQS>`.

We assume that the user has obtained a self-consistent solution of the
Kohn-Sham equations with Elk (a full tutorial can be found here :ref:`Elk SVO tutorial <SrVO3_elk>`). Also, the user needs to be familiar with the main in/output files of Elk, and how to run
the DFT code. Further information about Elk can be found on the :ref:`official Elk website <http://elk.sourceforge.net/>`.

Conversion for the DMFT self-consistency cycle
----------------------------------------------

Once the user has obtained the groundstate calculation, they will have to rerun Elk but with some small changes to the inputs in the elk.in file which will be explained below. The Elk part of the interface calculates and outputs the Wannier projectors. All downfolding related flags are set in the elk input file, and Elk determines automatically by symmetry equivalent sites. The TRIQS Elk converter then reads in these projectors along with the other Elk ascii files which would have been generated in the Elk ground state calculation. This information is then packed into the HDF5 file.

In the following, we use SrVO3 as an example to explain the flags required in the elk.in input file. An example elk.in of SrVO3 is available in the :ref:`SrVO3 tutorial <elk_SVO>`.

.. literalinclude:: ../tutorials/images_scripts/elk.in_SVO


The projectors are generated in Elk using these alterations in the elk.in file::
  
  task
   805

  wanproj
   1            !1) number of projector species
   2 2 3        !2) species number, l value and lm submatrix size
   7 8 9        !3) cubic harmonic lm indices 
   -0.055 0.07  !4) Correlated energy window in Hartrees (this is [-1.5, 1.9] in eV)

The first flag "task" specifies Elk to do the Wannier projector calculation in the Elk input convention. The "wanproj" flag specifies the information needed to generate the desired projectors (the exclamation marks are comment flags in Fortran). Below gives the meaning of each line:

#. The number of different species to generate the projectors for. If the material has multiple atoms of the same species, then the projectors will be generated for all of these atoms. Information about whether these atoms are symmetrically equivalent is written to the PROJ.OUT file. Generating projectors for multiple species requires lines 2) and 3) to be repeated for each projector species. 
#. The desired species index (relative to the order of the "atoms" input in elk.in), the l value and the number of wanted lm orbitals for the projectors. 
#. The lm indices (in the cubic harmonics by default) where lm = 1 refers to the s orbital, lm = 2,3,4 refers to the p orbitals, lm = 5,6,7,8,9 refers to the d orbitals and finally lm = 10,11,12,13,14,15,16 refers to the f orbitals. In the example above, this specifies that we want the t2g orbitals which has size 3 [last number in line 2)] and the lm indices are 7, 8 and 9 as specified in line 3). If, on the other hand, the user wishes to use all of the d-orbitals then all 5 orbital indices need to be included along with specifying that all 5 orbitals will be used [at the end of line 2)].
#. The correlated energy window (in Hartree) to generate the projectors within.

It should be noted that the indices in line 3) will change if another lm basis is used. 
The default is the cubic harmonic basis. The flags in elk.in required to change the spherical harmonic basis are::

  cubic
   .true.

  lmirep
   .true.

Above are the default inputs. If both of these flags are set to .false., the projectors will be generated in the complex spherical harmonic basis. It is possible to generate the projectors in Elk's irreducible basis by setting cubic to .false., but this is experimental and the TRIQS side of the interface is 
currently unable to convert the projectors in that basis. Finally, these projectors are written to file in the complex spherical harmonic basis.

Also, the input flag called "wanind", when set to .true. in elk.in, enables the user to input the lower and upper band indices (respectively) in line 4) of "wanproj" instead of the correlated energy window energies. An example of the inputs in the elk.in file::

  wanind
   .true.

  wanproj
   1      !No. of projector species
   2 2 3  !species, l, lm submatrix size
   7 8 9  !lm indices
   21 25  !correlated energy window band indices.

Note that for magnetic systems (spin uncoupled calculations), only the indices of the majority spin bands are used as inputs. The code calculates the indices for the minority spin automatically.


.. _Elk_files:

The rest of the elk.in file can remain unchanged. This 805 task calculates the projectors which are written into the WANPROJ_L**_S**_A****.OUT file(s). Note that the projectors will be written in the complex spherical harmonic basis. Along with this, the other written file (PROJ.OUT) specficies some information about the projectors (like atom equivalency, lm indices and so on) needed for reading the files into the TRIQS library. The PROJ.OUT file contains comments about its outputs. Here's a list of all of the input files needed for this part of the TRIQS converter:

#. WANPROJ_L**_S**_A****.OUT - file containing the projectors and band window indices.
#. PROJ.OUT -  specficies some information about the projectors.
   (like atom equivalency, lm indices and so on) needed for reading into the TRIQS library.
#. EIGVAL.OUT - contains the energies and latice vector coordinates for each k-point.
#. EFERMI.OUT - contains the Fermi energy.
#. KPOINTS.OUT - contains the k-point weights and lattice vectors.
#. SYMCRYS.OUT - has the crystal symmetries used for symmetries observables.
#. LATTICE.OUT - has lattice-Cartesian basis transformation matrices.
#. GEOMETRY.OUT - file with the lattice positions of every atom of each species.

Moreover, the Wannier charge density matrix (in WANCHARGE.OUT) and the Wannier spectral function (in WANSF_L**_S**_A****_0*.OUT) are calculated. These files are not used in the interface. 

As a side note, there are two other tasks which also generate the Wannier projectors. Task 804 generates the same outputs as 805 except it doesn't calculate the Wannier charge and spectral function. Task 806 outputs the same information as 805, but it generates the projectors (and other k-dependent variables) on a different user defined ngridk mesh. These tasks are parallelized with both OpenMP and MPI.

The Elk outputs are read into the TRIQS library using the following lines::

  from triqs_dft_tools.converters.elk import *
  Converter = ElkConverter(filename=filename, repacking=True)
  Converter.convert_dft_input()

The first two lines import and load the Elk converter module and the last line executes the conversion into the filename.h5 HDF5 file.


Data for post-processing - Correlated Spectral functions
--------------------------------------------------------

In case you want to do post-processing of your data using the module :class:`SumkDFTTools <dft.sumk_dft_tools.SumkDFTTools>`, some more files have to be converted to the HDF5 archive. Some of this has been laid out in :ref:`analysis`, but the Elk specific functions are described in the following sections. Below, we will discuss how to input the information for correlated spectral functions (band structures) and which is then calculated using "spaghettis" in :ref:`analysis`. However, band character band structure plots have not yet been implemented.

In Elk, the elk.in requires the altered flag::

  task
   820

As well as the wanproj flag (which has be discussed previously) and the plot1d flag (refer to the Elk manual). The new output files of this task required for the interface are:

#. BAND.OUT - the energy eigenvalues along the user specified path.
#. PROJ_WANBAND.OUT - same as PROJ.OUT but for band structure projectors
#. WANPROJ_L**_S**_A****_WANBAND.OUT - same as WANPROJ_L**_S**_A****.OUT
    but for band structure projectors.

(The ascii output files which have the new extension of _WANBAND.OUT are specific for this post processing calculation.)

The band structure information is converted into TRIQS by using::

  Converter.convert_bands_input()

Spectral function from Elk inputs
---------------------------------

Elk does not calculate the theta projectors for partial DOS calculations. Instead, Elk outputs the band characters into the file BC.OUT when using the elk.in task::

  task 
  803 

The contents of BC.OUT need to be converted into the HDF5 file by using the Elk Converter module::

  from triqs_dft_tools.converters.elk import *
  Converter = ElkConverter(filename=filename, repacking=True)
  Converter.dft_band_characters()

Once these have been saved to the HDF5 file (called "filename" here), the spectral function can be calculated with::

  SK.elk_dos(broadening=0.0, with_Sigma=True, with_dc=True, pdos=False, nk=None)

This outputs the total spectral function and the partial spectral function if enabled. Most of the user inputs are similar to the "SK.dos_parproj_basis()" module in :ref:`analysis`. The "pdos" flag when "True" enables the partial dos of each lm value to be calculated. It should be noted that these band characters are in Elk's irreducible lm basis and as such, the user has to check the irreducible representation used in Elk. This information can be found in the file ELMIREP.OUT after running task 10 (the DOS calculating task). The "nk" flag enables the calculation of the occupied spectral funciton. Here, nk needs to be the occupation density matrix (calculated from integrating the Green's function on the Matsubara axis) in the Bloch basis. This input needs to be in the same format as the occupation density matrix "deltaN" calculated in the sumk_DFT.calc_density_correction(dm_type='elk') module.


Spectral function Contour Plots (Fermi Surfaces) from Elk inputs
---------------------------------------------------------------

Here, we will discuss how to plot the Fermi surface contour or any other non-zero omega spectral function contour plot. This is currently tailored for the Elk inputs. From this point, we will refer to these contours as Fermi surfaces. The energy eigenvalues, projectors and so on required for the Fermi surface plot needs to be outputed from Elk. This is done by using::

  task 
  807 
  
in Elk, but unlike the previous Elk interface tasks, the k-mesh grid needs to be specified. This is done like using the same inputs as the Fermi surface calculations in Elk. In Elk, The user needs to specify the "plot3d" input flag used to generate the k-mesh which the interface variables are evaluated on. A simple example is for SrVO3 where plot3d would look something like::
  
  plot3d
  0.0 0.0 0.0 !1) origin
  1.0 0.0 0.0 !2) vertex 1
  0.0 1.0 0.0 !3) vertex 2
  0.0 0.0 1.0 !4) vertex 3
  32 32 32    !5) k-mesh grid size

Lines 1) to 4) specifies the corners (in lattice coordinates) of the k-grid box and line 5) is the grid size in each direction (see the Elk manual). If the user desires to plot a 2D plane, then the user should define the plane using lines 2) and 3) [relative to line 1)] and define line 4) to be the cross-product of lines 2) and 3) [i.e. the vector in line 4) is normal to the 2D plane]. The outputs will be in terms of the k-dependent quantities in the irreducible Brillouin zone (IBZ). The files needed for the interface are:

#. EIGVAL_FS.OUT - same as EIGVAL.OUT but the output is of the Fermi surface calculation.
#. KPOINT_FS.OUT - same as KPOINT.OUT but the output is of the Fermi surface calculation.   
#. PROJ_FS.OUT - same as PROJ.OUT but the output is of the Fermi surface calculation.   
#. WANPROJ_L**_S**_A****_FS.OUT - same as WANPROJ_L**_S**_A****.OUT but the output is of the Fermi surface calculation.   
#. EFERMI.OUT - contains the Fermi energy.
#. SYMCRYS.OUT - has the crystal symmetries used for symmetries observables.
#. LATTICE.OUT - has lattice-Cartesian basis transformation matrices.

(The ascii output files which have the extension of _FS.OUT are specific for this post processing calculation.)

These outputs are converted to the HDF5 file by::

  from triqs_dft_tools.converters.elk import *
  Converter = ElkConverter(filename=filename, repacking=True)
  Converter.convert_fs_input()

The spectral function for the Fermi surface plots are calculated with::

  SK.fs_plot(broadening=0.0, mesh=None, FS=True, plane=True, sym=True, orthvec=None, with_Sigma=True, with_dc=True)

The new flags specify the following: 

#. "FS" - determines whether the output will be the Fermi surface and uses the closest omega value to 0.0 in the mesh.  
#. "plane" - required to specify whether the Elk input parameters were generated on a k-mesh plane.
#. "sym" - needed if the IBZ will be folded out by using symmetry operations.
#. "orthvec" - (numpy array of length 3) needs to be specified if using "plane" as this input is the orthonormal vector to the 2D plane required for the folding out process.

To give the user a range of output capabilities, This routine can be used in the following ways:

#. If using "with_Sigma", the mesh will be the same as the self-energy. However, by setting FS=False, the user can input a mesh option if they desire the "Fermi surface" plots for each omega value (commensurate with the self-energy mesh) within the input range.

#. If the user is generating the DFT Spectral function plot (i.e. with_Sigma and with_dc both set to False), a mesh needs to be specified if FS=False. This function will output the spectral functions for the input mesh. Otherwise if FS is True, this would return the spectral function at omega=0.0.

The output files will have the form of "Akw_FS_X.dat" (X being either up, down or ud) if FS=True or "Akw_X_omega_Y.dat" (Y being the omega mesh index) otherwise. The latter file will have the omega values within the file (the fourth column). The first three columns of both output file types specifies the cartesian lattice vector (kx, ky, kz) and the last column is the spectral function values.


DFT+DMFT wavefunction dependent quantities
------------------------------------------

The DFT+DMFT wavefunctions and occupations are generated in Elk by diagonalizing the full DFT+DMFT density matrix at each k-point during task 808. These are used to update the electron density which is then used to solve the Kohn-Sham equations once to complete a FCSC DFT+DMFT cycle. It is possible to calculate quantities which are solely dependent on the wavefunctions and occupations by using the aforementioned diagonalized set. After generating the DMATDMFT.OUT file, these diagonalized wavefunctions and occupations are calculated by amending elk.in with the new task::

  task
   809

This will write the new second variational eigenvectors and occupations in the EVECSV.OUT and OCCSV.OUT binary files respectively. Then the wavefunction dependent quantities implemented within Elk can be calculated using these DFT+DMFT wavefunctions and occupations. Note that this only works for energy independent quantities. Also, the user has to ensure that the second variational eigenvectors will be used in determining the wavefunction dependent quantities. This is done by looking for the variable "tevecsv" in init0.f90 of Elk's source code and making sure that this is set to .true. for the task number the user wishes to use.



