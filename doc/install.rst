.. highlight:: bash

.. _install:


Packaged Versions of DFTTools
=============================

.. _ubuntu_debian:
Ubuntu Debian packages
----------------------

We provide a Debian package for the Ubuntu LTS Versions 16.04 (xenial) and 18.04 (bionic), which can be installed by following the steps outlined :ref:`here <triqslibs:triqs_debian>`, and the subsequent command::

        sudo apt-get install -y triqs_dft_tools

.. _anaconda:
Anaconda (experimental)
-----------------------

We provide Linux and OSX packages for the `Anaconda <https://www.anaconda.com/>`_ distribution. The packages are provided through the `conda-forge <https://conda-forge.org/>`_ repositories. After `installing conda <https://docs.conda.io/en/latest/miniconda.html>`_ you can install DFTTools with::

        conda install -c conda-forge triqs_dft_tools

See also `github.com/conda-forge/triqs_dft_tools-feedstock <https://github.com/conda-forge/triqs_dft_tools-feedstock/>`_.

.. _docker:
Docker
------

A Docker image including the latest version of DFTTools is available `here <https://hub.docker.com/r/flatironinstitute/triqs>`_. For more information, please see the page on :ref:`TRIQS Docker <triqslibs:triqs_docker>`.


Compiling DFTTools from source
==============================

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.
#. Likely, you will also need at least one impurity solver, e.g. the :ref:`CTHYB solver <triqscthyb:welcome>`.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/dft_tools`` repository from GitHub::

     $ git clone https://github.com/TRIQS/dft_tools dft_tools.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir dft_tools.build && cd dft_tools.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../dft_tools.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install


Installation steps for the use with WIEN2K version 14.2 and older
-----------------------------------------------------------------

.. warning::
   The following steps are only necessary if you use a Wien2k version
   14.2 and older. For newer versions (15.1 onwards), these things have
   been included already in the Wien2k package. There are no
   :program:`run_triqs` and :program:`runsp_triqs` scripts any more,
   since their functionality has been included in the Wien2k
   run scripts.

You need to take this last step manually since the Wien2k installation is not standard on all machines.
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

You will also need to insert manually a correct call of :file:`python` into
these scripts using an appropriate for your system MPI wrapper (mpirun,
mpprun, etc.), if needed.

Finally, you will have to change the calls to :program:`python_with_DMFT` to
your :program:`python` installation in the Wien2k :file:`path_to_Wien2k/run*` files.

 
Version compatibility
---------------------

Keep in mind that the version of ``dft_tools`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd dft_tools.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``dft_tools`` can be configured using CMake-options::

    cmake ../dft_tools.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_dft_tools|
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Check test coverage when testing                                | -DTEST_COVERAGE=ON                            |
| (run ``make coverage`` to show the results; requires the        |                                               |
| python ``coverage`` package)                                    |                                               |
+-----------------------------------------------------------------+-----------------------------------------------+
