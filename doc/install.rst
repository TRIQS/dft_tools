
.. highlight:: bash

Installation
============

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

#. Likely, you will also need at least one impurity solver, e.g. the :ref:`CTHYB solver <triqscthyb:welcome>`.

Installation steps 
------------------

#. Download the sources from github:: 
 
     $ git clone https://github.com/TRIQS/dft_tools.git src
 
#. Create an empty build directory where you will compile the code:: 
 
     $ mkdir build && cd build 
 
#. In the build directory call cmake specifying where the TRIQS library is installed:: 
 
     $ cmake -DTRIQS_PATH=path_to_triqs ../src 
 
#. Compile the code, run the tests and install the application:: 
 
     $ make 
     $ make test 
     $ make install 

Installation steps for use with WIEN2K
---------------------------------------

#. You need to take this last step manually since the Wien2k installation is not standard on all machines.
   After the above installation several files will be installed into::
  
     path_to_TRIQS_install_directory/share/triqs/Wien2k_SRC_files/SRC_templates
 
   These files are:

   * :file:`case.cf_f_mm2` and :file:`case.cf_p_cubic` containing matrices for
     the complex->cubic transformation of the local angular basis

   * :file:`case.indmftpr` is a template for the input file needed by the
     :program:`dmftproj` program. This program constructs a set of localized
     orbitals representing correlated states.

   These files then have to be copied manually to
   :file:`path_to_Wien2k/SRC_templates`, where :file:`path_to_Wien2k` is the path
   to the Wien2K main directory. 

   When building the Wien2k extension module, the :program:`dmftproj` is
   compiled and installed it into :file:`path_to_triqs/bin`. 

   In addition, :file:`path_to_Wien2k/SRC_templates` also contains
   :program:`run_triqs` and :program:`runsp_triqs` scripts for running Wien2k+DMFT
   fully self-consistent calculations. These files should be copied to
   :file:`path_to_Wien2k`, and set as executables by running::

     $ chmod +x run*_triqs 

   You will also need to insert manually a correct call of :file:`pytriqs` into
   these scripts using an appropriate for your system MPI wrapper (mpirun,
   mpprun, etc.), if needed. Search for *pytriqs* within the scripts to locate the
   appropriate place for inserting the :file:`pytriqs` call.

   Finally, you will have to change the calls to :program:`python_with_DMFT` to
   :program:`pytriqs` in the Wien2k :file:`path_to_Wien2k/run*` files.
 
Version compatibility 
--------------------- 
 
Be careful that the version of the TRIQS library and of the dft tools must be 
compatible (more information on the `TRIQS website 
<http://ipht.cea.fr/triqs/versions.html>`_). If you want to use a version of 
the dft tools that is not the latest one, go into the directory with the sources 
and look at all available versions:: 
 
     $ cd src && git tag 
 
Checkout the version of the code that you want:: 
 
     $ git co 1.0.0 
 
Then follow the steps 2 to 5 described above to compile the code. 
