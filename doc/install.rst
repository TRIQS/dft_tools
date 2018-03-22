.. highlight:: bash

.. _install:

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources from github::

     $ git clone https://github.com/triqs/app4triqs.git app4triqs.src

#. Create an empty build directory where you will compile the code::

     $ mkdir app4triqs.build && cd app4triqs.build

#. Make sure that you have added the TRIQS and Cpp2Py installation to your environment variables::

     $ source path_to_triqs/share/cpp2pyvarsh.sh
     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake::

     $ cmake ../app4triqs.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Be careful that the version of the TRIQS library and of the solver must be
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
As the TRILEX project is currently still in alpha phase (unstable), it can
only be compiled against the unstable branch of triqs.

Custom CMake options
--------------------

Functionality of ``app4triqs`` can be tweaked using extra compile-time options passed to CMake::

    cmake -DOPTION1=value1 -DOPTION2=value2 ... ../app4triqs.src

+-----------------------------------------------------------------------+-----------------------------------------------+
| Options                                                               | Syntax                                        |
+=======================================================================+===============================================+
| Specify an installation path other than path_to_triqs                 | -DCMAKE_INSTALL_PREFIX=path_to_app4triqs      |
+-----------------------------------------------------------------------+-----------------------------------------------+
| Build the documentation locally                                       | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------------+-----------------------------------------------+
