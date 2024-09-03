(changelog)=

# Changelog

## Version 3.3.1

DFTTools Version 3.3.1 is a patch release that restores compatibility
against recent numpy and scipy versions.

We thank all contributors: Sophie Beck, Alexander Hampel, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* fix compatibility against numpy 2.0
* replace numpy.lib.pad with numpy.pad
* use scipy.integrate.simpson instead of scipy.integrate.simps (#255)
* clean-up of optics part using wannierberri

### Doc
* update list of DFT codes on landing page


## Version 3.3.0

DFTTools Version 3.3.0 is a release that
* is compatible with TRIQS 3.3.x
* includes the latest app4triqs changes
* introduce `dc_imp_dyn` attribute in sumk object to store dynamic part of DC potential
* allows using MeshDLRImFreq as Sumk mesh
* improved standard behavior of block struct (#248) (see below for details)

We thank all contributors: Sophie Beck, Thomas Hahn, Alexander Hampel, Henri Menke, Dylan Simon, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* fix settings environment variables
* remove constrains on mpi size with vasp

### feat
* allow dict/np.ndarrays input in `symm_deg_gf`
* introduce `dc_imp_dyn` attribute in sumk object to store dynamic part of DC potential
* allows using MeshDLRImFreq as Sumk mesh
* previously the default `gf_struct_solver` in a initialized blockstructure had keys `up` / `down`, inconsistent with the default behavior after running `analyse_block_structure`: `up_0` / `down_0`. Now the default solver structure always has the `_0`
in the key.
* old behavior resulted in error when analyse_block_structure was called
twice
* to correctly use `analyse_block_structure` now use `extract_G_loc(transform_to_solver_blocks=False)`
* changed `density_matrix` function to use directly `extract_G_loc()` if `using_gf` is selected as option.
* print deprecation warning in `density_matrix`, since this can be achieved via `extract_G_loc` and `[G.density() for G in Gloc]`
* new function `density_matrix_using_point_integration()` for old point integration method option
* enforce in `analys_block_structure` that input dm or G is list with length of `n_corr_shells`
* correct doc string for how include_shells are given


### build
* bump actions/cache restore/save to version 4
* fix intel f2py build of elk converter (#249)
* fix MacOS X build: add ninja as req
* py312 add setuptools and meson as explicit dep for f2py
* add packaging directory to cmake and set version automatically

### doc
* add util module to autodoc


## Version 3.2.1

DFTTools Version 3.2.1 is a patch release that contains a few bug fixes. The following non breaking changes have been made:
* fix depracted scipy.compress depr -> numpy.compress
* fix incorrect numpy data type for Max OS ARM
* fix a bug in SumkDFT.calc_density_correction: see issue #250
* fix a bug in the Wannier90 Converter when the disentanglement window isn't set by the user (see issue #252)
* doc: fix typo in doi id of DC function

We thank all contributors: Sophie Beck, Alexander Hampel


## Version 3.2.0

DFTTools Version 3.2.0 is a release that
* is compatible with TRIQS 3.2.x
* introduces a coherent mesh of the SumK object passed on init used throughout all member functions (see below for details)
* unifies post-processing routines in `sumk_dft_tools` for all DFT codes to calculate DOS and spectral functions
* adds support for non-collinear projectors when using Vasp 6
* adds transport / optics calculation support with Elk
* introduces a new `dft_input` variable `dft_code` that identifies automatically the used dft code
* adds transport / optics calculation with wannier90 using WannierBerri
* moves all optics and transport related functions to a new module `sumk_dft_transport.py`
* new chemical potential finder (fully backward combatible), that supports besides bisection also gradient finders
* new double counting routines, moved to a separate function, improving spin-dependent formulas
* improve performance of all routines performing k-sums. `extract_G_loc` is up to 5 times faster now

We thank all contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Alyn James, Harrison LaBollita, Dario Fiore Mosca, Oleg Peil, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* general fixes to match new triqs 3.2 release
* SumK requires now to pass a mesh on init to clarify the mesh on which it operates
* rename / unify name of `sumk.Sigma_imp_iw` and `sumk.Sigma_imp_w` -> `sumk.Sigma_imp`
* remove `iw_or_w` arguments
* `sumk_dft_tools.py` rewritten to have single routines to calculate DOS (`dos_wannier_basis` renamed to `density_of_states`), spaghettis and (Elk specific for now) spectral contours
* occupied DOS can be calculated (`sumk_dft_tools.occupations()` is needed to be calculated first)
* analysis.rst and conv_elk.rst updated to improve routine descriptions and includes example figures
* remove any transport from `sumk_dft_tools.py` and move to `sumk_dft_transport`
* outsource `calc_DC_from_density` into util.py and cleanup
* Added new way to compute double counting. This is moved to a separate function and the old method is kept to modify in-place the double counting and keep compatibility. A legacy interface was kept for the old method (using integers to denote DC type). The new convention follow the notation introduced here (http://dx.doi.org/10.1038/s41598-018-27731-4)
* Added tests for DC calculation to compare with old implementation
* Added 2 new methods to find chemical potential, refactored DC calculation with stateless function while keeping legacy code
* In magnetic calculations, the dichotomy adjustment is struggling to find the mu (maxes out iterations). Added new methods to find the dft mu: newton (fastest but can fail) and brent (more stable)
* fix: fix f2py command for numpy ver >1.22
* fix: np.int / np.float / np. complex are  deprecated (np v1.20) / removed (np v1.24)
* Fixed DC formulas with SOC (#227)
* Fixes issue with projected A(k,w) (#225)
* fix obsolote iw_or_w param in calc_mu
* improve performance of extract Gloc by factor 5! huge cleanup
* edit SumKDFT class to take gf mesh at initialization and force Sigma to have that mesh
* update mpi.all_reduce calls 52bccac
* issue #216 correctly use beta when calling density on MeshReFreq

### Elk
* Elk Transport code and subsequent updates (#229)
* elk-interface bug fixes (#228)
* updated Elk tests and rewritten test scripts (.h5 files remain unchanged)
* New converter routines to read in Elk data for sumk_dft_tools.spectral_contours() (Elk k-mesh generator and checker needs to be optimized as it's currently slow). commented out Elk "bandcharacter" conversion from Elk converter and Elk DFT+DMFT PDOS code which used it (this method needs to be checked)

### Vasp
* change normion default to False
* change Vasp NiO tutorial scripts to reflect changes to sumk
* fix deprecated safeconfigparser
* fixes a bug in the Vasp charge self-consistent update step
* adjust NiO reference h5-file
* add test of NiO with two correlated shells
* fix tests of LuNiO3 and SrVO3 after changes
* fix mapping from shell/ions to corr-shells in converter

### w90 + QE
* feat: optical prop with Wannier90 and WannierBerri, see documentation for details
* fix bug for Gamma only mode
* w90 conv more generous matching to find fermi
* Updated W90 converter: bug fixes for SOC, code restructured, more tests for `add_lambda` and `bloch_basis`
* BUGFIX: changed character in QE output for reading occupations with
* Split wannier90 tests up
* Fix for wannier converter: reordering of orbital/spin order only necessary for vasp 5

### clean
* remove Gf indices and remove calc_dc_for_density (unused)


## Version 3.1.1

DFTTools Version 3.1.1 is a patch release that contains a few bug fixes.
In particular, we resolve incompatibility with recent numpy versions.

We thank all contributors: Alexander Hampel, Harry LaBollita

Find below an itemized list of changes in this release.

### General
* fix: np.int / np.float / np. complex are  depracted (np v1.20) / removed (np v1.24)
* Fixes issue with projected A(k,w) (#225)


## Version 3.1.0

DFTTools Version 3.1.0 is a release that
* is compatible with TRIQS 3.1.x
* includes a major update for the Wannier90 converter (see below for details)
* updates sumk_dft to allow for charge self-consistent DFT+DMFT calculations with Quantum Espresso (dm_type = 'qe')
* adds a indmftpr helper script to prepare the case.indmftpr file for the dmftproj program
* uses the latest [app4triqs/3.1.x](https://github.com/TRIQS/app4triqs) skeleton

### Wannier90 Converter
* allow for charge self-consistent DFT+DMFT calculations
* spin-orbit coupling implemented
* option to add a local spin-orbit term to t2g local Hamiltonian (for now just for a single impurity. Fixed in next version)
* additional choices and added checks for different bases (rot_mat): hloc_diag, wannier (already implemented previously), none
* code restructured, more tests
* MPI speedup of the Fourier transform
* added new test in w90_convert.py for rot_mat_type='hloc_diag'
* update documentation of W90 Converter
* bugfix: This fix makes the function find_rot_mat() safer to use in case there are errors in finding the correct mapping. The converter will now abort if the agreement in mapping is below a user-definable threshold.

### Change in gf_struct
* In line with TRIQS 3.1.x, the form of the Green's function's structure (`gf_struct`) has been modified (see [triqs changelog](https://triqs.github.io/triqs/latest/ChangeLog.html#change-in-gf-struct-objects) for more information)
* Instead of `gf_struct = [("up", [0, 1]), ("down", [0, 1])]`, the new convention uses `gf_struct = [("up", 2), ("down", 2)]`
* This modifies the form of `gf_struct_solver` (and `sumk`) in `block_structure` and `SumkDFT` as well.
* Backwards-compatibility with old, stored `block_structure` objects is given, however a warning is issued.
* A helper-function `triqs.gf.block_gf.fix_gf_struct_type(gf_struct_old)` is provided in triqs to manually bring `gf_struct`s to the new form.

### Documentation
* change to read the docs sphinx theme
* clean up various doc files
* use autosummary to build reference documentation
* update Vasp tutorials
* update Wannier90 documentation to reflect new features

### Cmake
* require triqs3.1+ in debian package dependencies
* bump required TRIQS Version to 3.1

### Other changes
* bugfix for analyse_block_structure in sumk_dft
* bugfix in blockstructure module for the case of #corr_shells != #ineq_shells
* fix float comparison tolerances and few minor things in tests
* Vasp Converter: fixed normalization of kwghts to allow symmetries
* bugfix in Elk converter when creating the symmetry matrices of low symmetry systems with multiple equivalent atoms
* vectorize various loops in dfttools
* fix various from_L_G_R calls that require now data layed out in C-order
* use nda over TRIQS_RUNTIME_ERROR in dos_tetra3d
* changed fermi weights from np array complex to float in accordance with h5 structure
* expose parameter max_loops in sum_k.calc_mu dichotomy

Thanks to all commit-contributors (in alphabetical order): Sophie Beck, Alexander Hampel, Alyn James, Jonathan Karp, Harry LaBollita, Max Merkel, H. L. Nourse, Hermann Schnait, Nils Wentzell,  @70akaline


## Version 3.0.0

DFTTools Version 3.0.0 is a major release that

* is compatible with TRIQS versions 3.0.x
* introduces compatibility with Python 3 (Python 2 is no longer supported)
* is now aligned with the general [app4triqs](https://github.com/TRIQS/app4triqs) application skeleton
* brings a major rework of the VASP interface, including thorough documentation, tutorials, a new Hamiltonian mode, the option to select bands instead of an energy window, and many small bugfixes.
* brings a major update of the block structure functionalities especially for SOC calculations, with detailed documentation and tutorials. Allows more control over the block structure coming from DFT, cutting out certain orbitals or throwing away off-diagonal elements when preparing input for the solver.
* New option in dmftproj to select the projection window using band indices instead of energie

### Restructuring

To be aligned with other applications for TRIQS, various files and folders had to be moved to new locations. The c++, fortran and python parts all are now in separate folders. The converter files have been more logically split into their own folders and name spaces. For example the Vasp converter is now located under `python/triqs_dft_tools/converters/vasp.py`. Especially the test folder structure was adapted to fit to the app4triqs skeleton, which separate folders for C++ and python tests.

### Dependency Management

We are managing the interdependencies of the various library components of triqs now using cmake.
Per default cmake will pull those dependencies from their corresponding
GitHub repositories, build them, and install these components together
with TRIQS, unless they are found in your system.
This behavior can be altered using the additional cmake options

* `-DBuild_Deps="Always"` - Always build dependencies, do not try to find them
* `-DBuild_Deps="Never"` - Never build dependencies, but find them instead

during the configuration step. See also the TRIQS documentation for more detailed instructions.

### Other Changes:

* Run port_to_triqs3 script
* Port py files to python3
* Update triqs python module name
* synchronize dfttools with app4triqs structure
* rename all h5 test archives according to test_name.ref.h5
* Changed 'orb' parameter to 'ish' for consistency in function summ_deg_gf in file sumk_dft.py
* small fix to read_inp_from_h5 function of Sumk
* renamed converters from app_converter.py to app.py
* look at the mesh of each shell of Sigma_imp, not just the first shell
* add function to find min and max of band energy, and add warning to set_Sigma if its mesh is smaller than the energy bounds
* warning for set_Sigma if ReFreqMesh is too small
* fixed a index bug that produced empty projectors for a unit cells with multiple shells in the VASP converter
* fixed a slicing bug for the calculation of the target density in the VASP converter, which selected 1 band less in the correlated window than required.
* added printout of complex part of local Hamiltonian in the Vasp converter
* doc on automatic basis rotations
* Bugfix in calculate_density_matrix for purely imaginary off-diagonals
* revamping the VASP interface documentation. Rewrote the interface with VASP guide. Removed the unused in doc/vasp. Start for SVO VASP tutorial as ipynb
* changed ref file for block structure test, since the order in dicts is not guaranteed the test failed as the order in py3 changed
* Vasp Converter: efermi is now read from LOCPROJ if DOSCAR does not contain it yet
* E-Fermi is read from DOSCAR not from LOCPROJ
* Added Tutorial for basis rotations: Sr2MgOsO6 w/o SOC
* Vasp converter add kpts and kpts_basis to h5
* many adjustments to Block structure and rotations including option to throw away certain parts of BlockGf
* implemented multiple ncsf VASP cycles
* Ignore imaginary part of the density when calculating mu
* Adjust hdf5 usage to changes in triqs
* Calculate diagonalization in solver blocks
* Do not use deprecated set_from_inverse_fourier
* add SOC tutorial
* add Block structure tutorial
* adding detailed Vasp tutorial
* Vasp converter now supports Hamiltonian mode
* Move setup files into separate bash scripts and adjust README
* Update README file with more detailed instructions


Thanks to all commit-contributors (in alphabetical order):
Markus Aichhorn, Alexander Hampel, Gernot Kraberger, Oleg Peil, Hermann Schnait, Malte Schueler, Nils Wentzell, Manuel Zingl


## Version 2.2.1

DFTTools Version 2.2.1 makes the application available
through the Anaconda package manager. We adjust
the install pages of the documentation accordingly.
We provide a more detailed description of the changes below.

* Add a section on the Anaconda package to the install page
* Add a LICENSE and AUTHORS file to the repository


## Version 2.2.0

* Ensure that the chemical potential calculations results in a real number
* Fix a bug in reading Wien2k optics files in SO/SP cases
* Some clarifications in the documentation
* Packaging/Jenkins/TRIQS/Installation adaptations

This is to a large extend a compatibility release against TRIQS version 2.2.0

Thanks to all commit-contributors (in alphabetical order):
Markus Aichhorn, Dylan Simon, Erik van Loon, Nils Wentzell, Manuel Zingl


## Version 2.1.x (changes since 1.4)

* Added Debian Packaging
* Compatibility changes for TRIQS 2.1.x
* Jenkins adjustments
* Add option to measure python test coverage
* VASP interface (and documentation)
* Added thermal conductivity in transport code (and documentation)
* BlockStructure class and new methods to analyze the block structure and degshells
* Multiple fixes of issues and bugs
* Major updates and restructuring of documentation



Thanks to all commit-contributors (in alphabetical order):
Markus Aichhorn, Gernot J. Kraberger, Olivier Parcollet, Oleg Peil, Hiroshi Shinaoka, Dylan Simon, Hugo U. R. Strand, Nils Wentzell, Manuel Zingl

Thanks to all user for reporting issues and suggesting improvements.
