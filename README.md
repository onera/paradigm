**ParaDiGM** (*Parallel Distributed General Mesh*) is a parallel computational geometry library under LGPL, with interfaces in C, Fortran and Python.

## Documentation  ##

User documentation is deployed on ONERA's internal GitLab pages server: https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev/index.html

## Build and install ##

### Basic Installation

Follow these steps to build **ParaDiGM** from the sources:

1. `git clone git@gitlab.onera.net:numerics/mesh/paradigm.git` (for ONERA users only)
2. `mkdir paradigm/build`
3. `cd paradigm/build`
4. `cmake ..`
5. `make`
6. `make install`
7. `./pdm_run` (if you want to run the test cases)

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

&#9888; *Note that shared libraries are necessary for Python interface (*`PDM_ENABLE_SHARED=ON`*).*

:warning:

#### Enable shared libraries
    PDM_ENABLE_SHARED=<ON | OFF> (default : ON)

#### Enable static libraries
    PDM_ENABLE_STATIC=<ON | OFF> (default : ON)

#### Enable ParMETIS library (parallel graph partitioning)
    PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find ParMETIS :

    PARMETIS_DIR=<path>

To link shared libraries, ParMETIS must be compiled with the `-fPIC` flag.<br>
CMake looks for
- `parmetis.h` and `metis.h` includes
- `parmetis` and `metis` libraries

#### Enable PT-Scotch library (parallel graph partitioning)

    PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find PT-Scotch :

    PTSCOTCH_DIR=<path>

To link shared libraries, PT-Scotch must be compiled with the `-fPIC` flag.<br>
CMake looks for
- `ptscotch.h` include file
- `scotch`, `scotcherr`, `ptscotch`, `ptscotcherr` libraries

#### Enable long global numbers
    PDM_ENABLE_LONG_G_NUM=<ON | OFF> (default : ON)
- `ON`  : `PDM_g_num_t` type is `long int`
- `OFF` : `PDM_g_num_t` type is `int`

#### Enable documentation compilation
    PDM_ENABLE_DOC=<ON | OFF> (default : OFF)
Once built, the documentation can be found in `build/doc/sphinx/html` and launch `index.html` file

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

## Issues ##

Issues can be reported directly on the [Issues](https://gitlab.onera.net/numerics/mesh/paradigm/-/issues) section.

## License ##

**ParaDiGM** is available under the LGPL3 license (https://www.gnu.org/licenses/lgpl-3.0.fr.html).
