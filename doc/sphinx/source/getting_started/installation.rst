.. _installation:

Installation
############

Dependencies
============

General dependencies for building **ParaDiGM** are:

  * a C++ compiler, Fortran is optional
  * `CMake <https://cmake.org/>`_ (version 3.16 or higher)
  * An **MPI** distribution (**required**).

Optional Dependencies
---------------------

Depending on the enabled options, the following dependencies may be required:

* **Python 3.8+** (with the **NumPy** module): Required if :code:`PDM_ENABLE_PYTHON_BINDINGS` is :code:`ON`.
* **Mpi4Py**: Required if :code:`PDM_ENABLE_PYTHON_BINDINGS` is :code:`ON`.
* **PT-Scotch**: Used for parallel graph partitioning.
* **ParMETIS**: Used for parallel graph partitioning.
* **Doxygen**: Required to generate C++/Fortran API documentation.
* **Breathe** and **Sphinx**: Required to compile the final documentation (if :code:`PDM_ENABLE_DOC` is :code:`ON`).
* **BLAS/LAPACK**: Used if :code:`PDM_ENABLE_BLASLAPACK` is :code:`ON`.
* **CUDA**: Required if :code:`PDM_ENABLE_CUDA` is :code:`ON`.


Basic Installation
==================

Follow these steps to build **ParaDiGM** from the sources:

.. code-block:: sh

  mkdir build
  cd build
  cmake ..
  make
  make install


If installation fails, use the following CMake options.


CMake general options
=====================

CMake options are passed on the command line in the format:

.. code-block:: sh

    cmake -D<option1_name>=<option1_value> ... -D<option2_name>=<option2_value>

**Build Type**

Defines the optimization and debugging level.

.. code-block:: sh

    CMAKE_BUILD_TYPE=<TYPE> (default: Release)

Common options are: :code:`Debug`, :code:`Release`, :code:`RelWithDebInfo`, :code:`MinSizeRel`, :code:`Sanitize`.

**Installation prefix**

Defines the path where the compiled files will be installed by :code:`make install`.

.. code-block:: sh

  CMAKE_INSTALL_PREFIX=<prefix>


Paradigm Specific Options
=========================

The following options control the features and library build types.

Library Build Options
---------------------

**Build Shared Libraries**

.. code-block:: sh

    PDM_ENABLE_SHARED=<ON | OFF> (default: ON)

.. note::
    If :code:`PDM_ENABLE_PYTHON_BINDINGS` is :code:`ON`, this option is **forced to ON**.

**Build Static Libraries**

.. code-block:: sh

    PDM_ENABLE_STATIC=<ON | OFF> (default: OFF)

**Enable Long Global IDs**

Controls the integer size used for global numbers (element or node IDs).

.. code-block:: sh

    PDM_ENABLE_LONG_G_NUM=<ON | OFF> (default: ON)

* :code:`ON` : The type :code:`PDM_g_num_t` is :code:`long int`.
* :code:`OFF` : The type :code:`PDM_g_num_t` is :code:`int`.

Interfaces and Extensions
-------------------------

.. _enable_fortran_interface:

**Enable Fortran Interface**

.. code-block:: sh

    PDM_ENABLE_Fortran=<ON | OFF> (default: ON)

**Build with Fortran MPI Module**

Enables Fortran preprocessing for using the standard Fortran MPI module (requires :code:`PDM_ENABLE_Fortran` to be :code:`ON`).

.. code-block:: sh

    PDM_ENABLE_Fortran_MPI_MODULE=<ON | OFF> (default: OFF)

.. _enable_python_interface:

**Enable Python Bindings**

Enables the compilation of the Python module.

.. code-block:: sh

    PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default: OFF)

If autodetection fails for **Python** and **NumPy**, use the variables from `FindPython` (see CMake documentation). For **Mpi4Py**, the :code:`Mpi4Py_DIR` variable may be helpful.

**Enable PDMA Extension (Advanced functionalities)**

.. code-block:: sh

    PDM_ENABLE_EXTENSION_PDMA=<ON | OFF> (default: OFF)

**Enable Anisotropic Agglomeration (if PDMA is enabled)**

This option is only available if :code:`PDM_ENABLE_EXTENSION_PDMA` is :code:`ON`.

.. code-block:: sh

    PDM_ENABLE_ANISO_AGGLO=<ON | OFF> (default: ON)

External Dependencies and Paths
-------------------------------

**Enable the use of BLAS / LAPACK (Linear algebra)**

.. code-block:: sh

    PDM_ENABLE_BLASLAPACK=<ON | OFF> (default: OFF)

If autodetection fails, you can use the variables from `FindBLAS` and `FindLAPACK`.


.. _parmetis: https://github.com/KarypisLab/ParMETIS
.. |parmetis| replace:: **ParMETIS**

**Enable the use of** |parmetis|_ **(Parallel graph partitioning)**

ParMETIS is searched for automatically. If found, the modern targets **ParMETIS::ParMETIS** and **Metis::Metis** are created and linked.

If autodetection fails, you must manually specify the root installation directory. The custom search module checks the following **CMake variables** (passed via :code:`-D`) and **environment variables** (exported in the shell) to locate the installation path (the directory containing :code:`include/` and :code:`lib/`):

.. code-block:: sh

    # CMake Variables (highest priority)
    PARMETIS_ROOT=<path>
    PARMETIS_DIR=<path>

    # Environment Variables (exported in the shell before calling CMake)
    # Example: export PARMETIS_DIR=/path/to/parmetis
    $ENV{PARMETIS_DIR}
    $ENV{PARMETIS_ROOT}

The search requires both the **ParMETIS** and **METIS** libraries to be found. If METIS is not located alongside ParMETIS, you can specify its path separately using :code:`METIS_ROOT` or :code:`METIS_DIR`.


.. _ptscotch: https://gitlab.inria.fr/scotch/scotch
.. |ptscotch| replace:: **PT-Scotch**

**Enable the use of** |ptscotch|_ **(Parallel graph partitioning)**

PT-Scotch is searched for automatically. If found, the primary target **PTScotch::PTScotch** and the sequential target **Scotch::Scotch** are linked.

If autodetection fails, you must specify the root installation directory. The custom search module checks the following **CMake variables** (passed via :code:`-D`) and **environment variables** (exported in the shell) to locate the installation path:

.. code-block:: sh

    # CMake Variables (highest priority)
    PTSCOTCH_ROOT=<path>
    PTSCOTCH_DIR=<path>

    # Environment Variables (exported in the shell before calling CMake)
    $ENV{PTSCOTCH_DIR}
    $ENV{PTSCOTCH_ROOT}
    $ENV{SCOTCH_DIR} (as a fallback)
    $ENV{SCOTCH_ROOT} (as a fallback)

The search requires the **ptscotch.h** include file and the main **ptscotch** library.


**Build with CUDA Support**

.. code-block:: sh

    PDM_ENABLE_CUDA=<ON | OFF> (default: OFF)

### Testing and Documentation

**Enable Sphinx Documentation Compilation**

Requires **Doxygen**, the **Breathe** Sphinx extension, and **Sphinx** installed.

.. code-block:: sh

    PDM_ENABLE_DOC=<ON | OFF> (default: OFF)

Once built, the documentation can be found in :code:`build/doc/sphinx/html/index.html`.

.. warning::
    Documentation **cannot** be built if :code:`CMAKE_BUILD_TYPE` is set to :code:`Sanitize`.

**Enable CTest Execution**

.. code-block:: sh

    PDM_ENABLE_TESTS=<ON | OFF> (default: ON)

**Enable Unit Test Compilation**

Requires the **doctest** framework (fetched via `FetchContent`).

.. code-block:: sh

    PDM_ENABLE_UNIT_TEST=<ON | OFF> (default: OFF)

**Enable Training Jupyter Notebooks**

.. code-block:: sh

    PDM_ENABLE_TRAINING=<ON | OFF> (default: OFF)


Compiler Selection
==================

You can specify compilers using environment variables when calling CMake:

.. code-block:: sh

    CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

Or by using the following CMake options:

.. code-block:: sh

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>


MPI CMake Options
=================

The project **requires** MPI and uses the standard CMake module **FindMPI.cmake**. Automatic detection based on environment and common paths is performed.

If automatic detection fails (e.g., if you are not using standard MPI wrappers like :code:`mpicc`), you must manually provide the paths to the MPI wrappers, which is the standard fallback mechanism for `FindMPI`:

**Specify MPI Wrappers (high-level FindMPI mechanism)**

.. code-block:: sh

    MPI_C_COMPILER=<C MPI wrapper>
    MPI_CXX_COMPILER=<CXX MPI wrapper>
    MPI_Fortran_COMPILER=<Fortran MPI wrapper>

Refer to the official `FindMPI documentation <https://cmake.org/cmake/help/latest/module/FindMPI.html>`_ for a complete list of variables and a detailed explanation of the detection process.
