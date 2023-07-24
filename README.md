**ParaDiGM** (*Parallel Distributed General Mesh*) is a parallel computational geometry library under LGPL.

## Documentation  ##

User documentation is deployed on the Gitlab pages server: https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev/index.html

## Build and install ##

### Basic Installation

Follow these steps to build ParaDiGM from the sources:

1. `git clone git@gitlab.onera.net:numerics/mesh/paradigm.git` (for ONERA users)
2. `mkdir paradigm/build`
3. `cd paradigm/build`
4. `cmake ..`
5. `make`
6. `make install`
7. `./pdm_run` (if you want to run the test cases)

### CMake general options

cmake -D<option1_name>=<option1_value> ... -D<option_name>=<option_value>

Prefix :
    CMAKE_INSTALL_PREFIX=<prefix>

Enable Fortran interface :
    PDM_ENABLE_Fortran=<ON | OFF> (default : OFF)

Enable Python interface :
    PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)
      If a simple autodetection fails, you can use these options to find Python :
        Python_ROOT_DIR=<path>
        Python_LIBRARY=<path>
        Python_INCLUDE_DIR=<path>
        Python_EXECUTABLE=<path>

      Refere to FindPython in the CMake documentation for more information.
      shared libraries are necessary for Python interface (PDM_ENABLE_SHARED=ON)

Enable shared libraries :
    PDM_ENABLE_SHARED=<ON | OFF> (default : ON)

Enable static libraries :
    PDM_ENABLE_STATIC=<ON | OFF> (default : ON)

Enable ParMETIS library (parallel graph partition) :
    PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)
      If a simple autodetection fails, you can use these options to find ParMETIS :
        PARMETIS_DIR=<path>

     To link shared libraries, ParMETIS has to be compiled with "-fPIC" option.

     CMake looks for :
        - parmetis.h and metis.h includes
        - parmetis and metis libraries

Enable PTSCOTCH library (parallel graph partition) :
    PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)
      If a simple autodetection fails, you can use these options to find ParMETIS :
        PARMETIS_DIR=<path>

     To link shared libraries, PTSCOTCH has to be compiled with "-fPIC" option.

     CMake looks for :
        - ptscotch.h include file
        - scotch, scotcherr, ptscotch, ptscotcherr libraries

Enable long global numbers
     PDM_ENABLE_LONG_G_NUM=<ON | OFF> (default : ON)
       - ON  : PDM_g_num_t type is "long int"
       - OFF : PDM_g_num_t type is "int"

Enable documentation mode :
     PDM_ENABLE_DOC=<ON | OFF> (default : OFF)
     Once build, the documentation can be found in build/doc/sphinx/html and launch `index.html` file

### Compiler choice

CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>

### CMake MPI options

    MPI_C_COMPILER=<C mpi wrapper>
    MPI_CXX_COMPILER=<CXX mpi wrapper>
    MPI_Fortran_COMPILER=<Fortran mpi wrapper>

If a simple autodetection fails, you can use these options to find MPI :

    MPI_<language>_LIBRARIES
    MPI_<language>_INCLUDE_PATH

Refere to FindMPI in the CMake documentation for more informations.

## Issues ##

Issues can be reported directly on [the Issues section](https://gitlab.onera.net/numerics/mesh/paradigm/-/issues).

## License ##

`ParaDiGM` is available under the LGPLv3 license (https://www.gnu.org/licenses/lgpl-3.0.fr.html).
