# ------------------------------------------------------------------------------
# General CMake settings
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
cmake_policy(SET CMP0074 NEW) # force find_package to take <PackageName>_ROOT variables into account

if(NOT DEFINED PROJECT_ROOT)
  set(PROJECT_ROOT ${CMAKE_SOURCE_DIR} CACHE PATH "Root directory, where the submodules are populated")
endif()

#------------------------------------------------------------------------------
# Internal CMake module (shared with ParaDiGMA or CWIPI)
#------------------------------------------------------------------------------
set(PDM_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${PDM_CMAKE_DIR}/modules")

# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/version.cmake)
include(${CMAKE_MODULE_PATH}/find_subdir.cmake)
project(PDM VERSION ${PDM_DEF_VERSION_MAJOR} LANGUAGES Fortran C CXX)

# A recuperer depuis git ? --> Fail on recup la 2.1 --> A discuter avec Romain/Julien/Bérenger

set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

## Compiler flags
### C++ standard
set(SUPPORTED_CXX_STANDARDS 11 14 17 20)


# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
## System dependencies ###
if (NOT TARGET Python::Python OR NOT TARGET Python::NumPy)
  find_package(Python 3.8 REQUIRED COMPONENTS Development NumPy)
endif()
if (NOT TARGET MPI::MPI_CXX)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if (NOT TARGET Mpi4Py::Mpi4Py)
  find_package(Mpi4Py REQUIRED)
endif()

find_package(PTScotch)
find_package(ParMETIS)
find_package(BLAS)
find_package(LAPACK)
find_package(Doxygen)


# ------------------------------------------------------------------------------
# Features
# ------------------------------------------------------------------------------

option(PDM_ENABLE_LONG_G_NUM      "Enable long for global number" ON )
option(PDM_ENABLE_SHARED          "Build Shared Libraries"        ON )
option(PDM_ENABLE_STATIC          "Build Static Libraries"        OFF)

option(PDM_ENABLE_Fortran         "Build Fortran Interface"       ON )
option(PDM_ENABLE_PYTHON_BINDINGS "Build Python module"           OFF)

option(PDM_ENABLE_UNIT_TEST       "Enable unit test "             OFF)

if (PDM_ENABLE_Fortran_MPI_MODULE)
  set(CMAKE_Fortran_PREPROCESS ON)
  set(PDM_HAVE_FORTRAN_MPI_MODULE 1)
endif()

if (PDM_ENABLE_LONG_G_NUM)
  set(PDM_LONG_G_NUM 1)
endif ()

configure_file(pdm_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")

# ------------------------------------------------------------------------------
# Compilation
# ------------------------------------------------------------------------------
set(src_dir     ${PROJECT_SOURCE_DIR}/src)
set(test_dir    ${PROJECT_SOURCE_DIR}/test)
set(include_dir ${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE all_src_files
  CONFIGURE_DEPENDS "${src_dir}/*.cxx" "${src_dir}/*.cpp" "${src_dir}/*.c"
)

file(GLOB_RECURSE all_src_fortran_files
  CONFIGURE_DEPENDS "${src_dir}/*.F90"
)

get_all_subdirectories(PDM_SUBDIRS ${src_dir})

set(src_files ${all_src_files})
list(FILTER src_files EXCLUDE REGEX ".*\\.test\\.cpp$")

set(LIB_TYPES "")
if(PDM_ENABLE_STATIC)
  list(APPEND LIB_TYPES STATIC)
endif()
if(PDM_ENABLE_SHARED)
  list(APPEND LIB_TYPES SHARED)
endif()

foreach(lib_type IN LISTS LIB_TYPES)

  add_library(pdm_${lib_type} ${lib_type} ${src_files})
  add_library(pdmf_${lib_type} ${lib_type} ${all_src_fortran_files})

  set_target_properties(pdm_${lib_type} PROPERTIES OUTPUT_NAME "pdm")
  set_target_properties(pdmf_${lib_type} PROPERTIES OUTPUT_NAME "pdmf")

  target_include_directories(pdm_${lib_type} PUBLIC
    ${PDM_SUBDIRS}
    $<BUILD_INTERFACE:${src_dir}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )

  # > Les 2 marches
  # foreach(subdir IN LISTS PDM_SUBDIRS)
  #   target_include_directories(pdm_${lib_type} PUBLIC
  #       $<BUILD_INTERFACE:${subdir}>
  #       $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  #   )
  # endforeach()

  target_link_libraries(pdm_${lib_type}
    PUBLIC
      MPI::MPI_CXX
  )

  target_compile_definitions(pdm_${lib_type} PRIVATE PDM_USE_MULTIPART=0)

  if (PTSCOTCH_FOUND)
    target_compile_definitions(pdm_${lib_type} PRIVATE PDM_HAVE_PTSCOTCH=1)
    target_link_libraries(pdm_${lib_type}
    PUBLIC
      PTScotch::PTScotch
      Scotch::Scotch
      Scotch::ScotchErr
    )
  endif()

  if (PARMETIS_FOUND)
    target_compile_definitions(pdm_${lib_type} PRIVATE PDM_HAVE_PARMETIS=1)
    target_link_libraries(pdm_${lib_type}
    PUBLIC
      ParMETIS::ParMETIS
    )
  endif()

  if (BLAS_FOUND)
    target_compile_definitions(pdm_${lib_type} PRIVATE PDM_HAVE_BLAS=1)
    target_link_libraries(pdm_${lib_type}
    PUBLIC
      BLAS::BLAS
    )
  endif()

  if (LAPACK_FOUND)
    target_compile_definitions(pdm_${lib_type} PRIVATE PDM_HAVE_LAPACK=1)
    target_link_libraries(pdm_${lib_type}
    PUBLIC
      LAPACK::LAPACK
    )
  endif()

  target_link_libraries(pdmf_${lib_type}
    PUBLIC
      pdm_${lib_type}
  )
endforeach()

# ----------------------------------------------------------------------
# Unit test
# ----------------------------------------------------------------------
if(PDM_ENABLE_UNIT_TEST)
  include(FetchContent)
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
  )
  FetchContent_GetProperties(doctest)
  if(NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_subdirectory(${doctest_SOURCE_DIR} ${doctest_BINARY_DIR})
  endif()
endif()

# ----------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------
include(CTest)
if (${PROJECT_NAME}_ENABLE_TESTS)
  add_subdirectory(test)
endif()

# ----------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()
# ----------------------------------------------------------------------
