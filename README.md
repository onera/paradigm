<p align="center">
    <img src="doc/_static/logo.svg" alt="Logo" width="25%"/>
</p>

**ParaDiGM** (*Parallel Distributed General Mesh*) is a library for managing meshes in a massively parallel distributed environment, with interfaces in C, Python and Fortran.

## 📖 Documentation

The user documentation is available either on [the GitHub page](https://onera.github.io/paradigm/index.html) or on [ONERA's internal Gitlab](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev/index.html).

## 🛠️ Build and install

### Dependencies

General dependencies for building **ParaDiGM** are:
- a C compiler
- [CMake](https://cmake.org/) (version 3.16 or higher)
- an MPI distribution

### Basic Installation

Follow these steps to build **ParaDiGM** from the sources:

First, clone the repository
  - either from [GitHub](https://github.com/onera/paradigm): `git clone git@github.com:onera/paradigm.git`
  - or from [GitLab](https://gitlab.onera.net/numerics/mesh/paradigm) (for ONERA users only): `git clone git@gitlab.onera.net:numerics/mesh/paradigm.git`

Then, use the following commands:
1. `cd paradigm`
2. `git submodule update --init extensions/paradigma` (if you want to enable [**ParaDiGMA**](https://gitlab.onera.net/numerics/mesh/paradigma))
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`
7. `make install`
8. `./pdm_run` (if you want to run the test cases)

<details>

<summary>Configuration with CMake</summary>

### CMake general options
    cmake -D<option1_name>=<option1_value> ... -D<option2_name>=<option2_value>

#### Installation prefix
    CMAKE_INSTALL_PREFIX=<prefix>

#### Enable Fortran interface
    PDM_ENABLE_Fortran=<ON | OFF> (default : OFF)

#### Enable Python interface
    PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)

If a simple autodetection fails, you can use these options to find Python :

    PYTHON_LIBRARY=<path>
    PYTHON_INCLUDE_DIR=<path>

Refer to [FindPython](https://cmake.org/cmake/help/latest/module/FindPython.html) in the CMake documentation for more information.

:warning: *Note that a shared library is necessary for Python interface:* `PDM_ENABLE_SHARED=ON`



#### Build shared library
    PDM_ENABLE_SHARED=<ON | OFF> (default : ON)

#### Build static library
    PDM_ENABLE_STATIC=<ON | OFF> (default : ON)


#### Build tests
    PDM_ENABLE_TESTS=<ON | OFF> (default : ON)

#### Enable the use of [ParMETIS](https://github.com/KarypisLab/ParMETIS) (parallel graph partitioning)
    PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find ParMETIS :

    PARMETIS_DIR=<path>

To link shared libraries, ParMETIS must be compiled with the `-fPIC` flag.<br>
CMake looks for
- `parmetis.h` and `metis.h` includes
- `parmetis` and `metis` libraries

#### Enable the use of [PT-Scotch](https://gitlab.inria.fr/scotch/scotch) (parallel graph partitioning)

    PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find PT-Scotch :

    PTSCOTCH_DIR=<path>

To link shared libraries, PT-Scotch must be compiled with the `-fPIC` flag.<br>
CMake looks for
- `ptscotch.h` include file
- `scotch`, `scotcherr`, `ptscotch`, `ptscotcherr` libraries

#### Enable the use of [BLAS](https://www.netlib.org/blas/) / [LAPACK](https://www.netlib.org/lapack/) (linear algebra)
    PDM_ENABLE_BLASLAPACK=<ON | OFF> (default : OFF)

#### Enable ``long int`` for global IDs
    PDM_ENABLE_LONG_G_NUM=<ON | OFF> (default : ON)
- `ON`  : `PDM_g_num_t` type is `long int`
- `OFF` : `PDM_g_num_t` type is `int`

#### Enable documentation compilation
    PDM_ENABLE_DOC=<ON | OFF> (default : OFF)
Once built (using the ``make doc`` command), the documentation home page can be found in `build/doc/sphinx/html/index.html`

### Compiler choice
    CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>

### CMake MPI options
    MPI_C_COMPILER=<C MPI wrapper>
    MPI_CXX_COMPILER=<CXX MPI wrapper>
    MPI_Fortran_COMPILER=<Fortran MPI wrapper>

If a simple autodetection fails, you can use these options to find MPI :

    MPI_<language>_LIBRARIES
    MPI_<language>_INCLUDE_PATH

Refer to [FindMPI](https://cmake.org/cmake/help/latest/module/FindMPI.html) in the CMake documentation for more information.

</details>

## 🚀 Quick start

Examples of use can be found in the ``test`` and ``training`` directories and in the [documentation](#documentation).


## 🤕 Issues

Issues can be reported directly in the [Issues](https://gitlab.onera.net/numerics/mesh/paradigm/-/issues) section.

## ⚖️ License

**ParaDiGM** is available under the [LGPL3 license](https://www.gnu.org/licenses/lgpl-3.0.en.html).

## © Copyright

Copyright 2026, ONERA The French Aerospace Lab