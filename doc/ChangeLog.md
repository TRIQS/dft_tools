Version 3.0.0
=============

DFTTools Version 3.0.0 is a major release that

* is compatible with TRIQS versions 3.0.x
* introduces compatibility with Python 3 (Python 2 is no longer supported)
* is now aligned with the general [app4triqs](https://github.com/TRIQS/app4triqs) application skeleton
* brings a major rework of the VASP interface, including thorough documentation, tutorials, a new Hamiltonian mode, the option to select bands instead of an energy window, and many small bugfixes.
* brings a major update of the block structure functionalities especially for SOC calculations, with detailed documentation and tutorials. Allows more control over the block structure coming from DFT, cutting out certain orbitals or throwing away off-diagonal elements when preparing input for the solver.
* New option in dmftproj to select the projection window using band indices instead of energie

Restructuring
-------------
To be aligned with other applications for TRIQS, various files and folders had to be moved to new locations. The c++, fortran and python parts all are now in separate folders. The converter files have been more logically split into their own folders and name spaces. For example the Vasp converter is now located under `python/triqs_dft_tools/converters/vasp.py`. Especially the test folder structure was adapted to fit to the app4triqs skeleton, which separate folders for C++ and python tests.

Dependency Management
--------------------
We are managing the interdependencies of the various library components of triqs now using cmake.
Per default cmake will pull those dependencies from their corresponding
GitHub repositories, build them, and install these components together
with TRIQS, unless they are found in your system.
This behavior can be altered using the additional cmake options

* `-DBuild_Deps="Always"` - Always build dependencies, do not try to find them
* `-DBuild_Deps="Never"` - Never build dependencies, but find them instead

during the configuration step. See also the TRIQS documentation for more detailed instructions.

Other Changes:
-------------
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


Version 2.2.1
=============

DFTTools Version 2.2.1 makes the application available
through the Anaconda package manager. We adjust
the install pages of the documentation accordingly.
We provide a more detailed description of the changes below.

* Add a section on the Anaconda package to the install page
* Add a LICENSE and AUTHORS file to the repository


Version 2.2.0
=============

* Ensure that the chemical potential calculations results in a real number
* Fix a bug in reading Wien2k optics files in SO/SP cases
* Some clarifications in the documentation
* Packaging/Jenkins/TRIQS/Installation adaptations

This is to a large extend a compatibility release against TRIQS version 2.2.0

Thanks to all commit-contributors (in alphabetical order):
Markus Aichhorn, Dylan Simon, Erik van Loon, Nils Wentzell, Manuel Zingl


Version 2.1.x (changes since 1.4)
=================================

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
