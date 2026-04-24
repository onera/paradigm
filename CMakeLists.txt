# ------------------------------------------------------------------------------
# General CMake settings
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
cmake_policy(SET CMP0074 NEW) # force find_package to take <PackageName>_ROOT variables into account

if (NOT DEFINED PROJECT_ROOT)
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

set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

## C++ standard
set(SUPPORTED_CXX_STANDARDS 11 14 17 20)

#------------------------------------------------------------------------------
# Default compiler flags
#------------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(${CMAKE_MODULE_PATH}/default_flags.cmake)
endif()

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
if (NOT TARGET MPI::MPI_CXX)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if (PDM_ENABLE_PYTHON_BINDINGS)
  if (NOT TARGET Python::Python OR NOT TARGET Python::NumPy)
    find_package(Python 3.8 REQUIRED COMPONENTS Development NumPy)
  endif()
  if (NOT TARGET Mpi4Py::Mpi4Py)
    find_package(Mpi4Py REQUIRED)
  endif()
endif()

find_package(PTScotch)
find_package(ParMETIS)
find_package(Doxygen)

# ------------------------------------------------------------------------------
# Features
# ------------------------------------------------------------------------------
option(PDM_ENABLE_EXTENSION_PDMA     "Enable advanced functionalities extension" OFF)
option(PDM_ENABLE_LONG_G_NUM         "Enable long for global number"             ON )
option(PDM_ENABLE_SHARED             "Build Shared Libraries"                    ON )
option(PDM_ENABLE_STATIC             "Build Static Libraries"                    OFF)
if (PDM_ENABLE_EXTENSION_PDMA)
  option(PDM_ENABLE_ANISO_AGGLO      "Enable anisotropic agglomeration"          ON )
endif()

option(PDM_ENABLE_Fortran            "Build Fortran Interface"                   ON )
option(PDM_ENABLE_Fortran_MPI_MODULE "Build with Fortran MPI module"             OFF)
option(PDM_ENABLE_PYTHON_BINDINGS    "Build Python module"                       OFF)
option(PDM_ENABLE_CUDA               "Build with CUDA"                           OFF)
option(PDM_ENABLE_BLASLAPACK         "Build with BLAS/LAPACK"                    OFF)
option(PDM_ENABLE_PEDANTIC           "Force pedantic compilation (for dev only)" OFF)

option(PDM_ENABLE_TESTS              "Make CTest run the tests"                  ON )
option(PDM_ENABLE_UNIT_TEST          "Enable unit test "                         OFF)
option(PDM_ENABLE_DOC                "Enable Sphinx documentation"               OFF)
option(PDM_ENABLE_TRAINING           "Enable training Jupyter Notebooks"         OFF)

if (PDM_ENABLE_LONG_G_NUM)
  set(PDM_LONG_G_NUM 1)
endif()

if (PDM_ENABLE_EXTENSION_PDMA AND PDM_ENABLE_ANISO_AGGLO)
  set (PDM_HAVE_ANISO_AGGLO ON)
endif()

if (PDM_ENABLE_Fortran_MPI_MODULE)
  set(CMAKE_Fortran_PREPROCESS ON)
  set(PDM_HAVE_FORTRAN_MPI_MODULE 1)
endif()

if (PDM_ENABLE_PYTHON_BINDINGS)
  if (NOT PDM_ENABLE_SHARED)
    set (PDM_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE)
  endif()
  find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Mpi4Py REQUIRED)
endif()

if (PDM_ENABLE_CUDA)
  enable_language(CUDA)
endif()

if (PDM_ENABLE_BLASLAPACK)
  find_package(BLAS)
  find_package(LAPACK)
endif()

if (PTSCOTCH_FOUND)
  set (PDM_HAVE_PTSCOTCH ON)
endif()
if (PARMETIS_FOUND)
  set (PDM_HAVE_PARMETIS ON)
endif()
if (BLAS_FOUND)
  set (PDM_HAVE_BLAS ON)
endif()
if (LAPACK_FOUND)
  set (PDM_HAVE_LAPACK ON)
endif()

set (install_configs "")
configure_file(pdm_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")
list (APPEND install_configs ${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h)
if (PDM_ENABLE_Fortran)
  if ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h" IS_NEWER_THAN "${CMAKE_CURRENT_BINARY_DIR}/pdm_configf.h")
    include (FortranizeFile)
    fortranize ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")
    list (APPEND install_configs ${CMAKE_CURRENT_BINARY_DIR}/pdm_configf.h)
  endif()
endif()

#------------------------------------------------------------------------------
# Fake dlclose path
#------------------------------------------------------------------------------
set(PDM_BINARY_DLCLOSE_DIR ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose)
add_subdirectory(script)

# ------------------------------------------------------------------------------
# Compilation
# ------------------------------------------------------------------------------
set(src_dir ${PROJECT_SOURCE_DIR}/src)
if (PDM_ENABLE_EXTENSION_PDMA)
  set(pdma_dir     ${PROJECT_SOURCE_DIR}/extensions/paradigma)
  set(pdma_src_dir ${pdma_dir}/src)
endif()

file(GLOB_RECURSE all_src_files
  CONFIGURE_DEPENDS "${src_dir}/*.cxx" "${src_dir}/*.cpp" "${src_dir}/*.c"
)

file(GLOB_RECURSE all_src_fortran_files
  CONFIGURE_DEPENDS "${src_dir}/*.F90"
)

get_all_subdirectories(PDM_SUBDIRS ${src_dir})

if (PDM_ENABLE_EXTENSION_PDMA)
  file(GLOB_RECURSE pdma_sources
    CONFIGURE_DEPENDS "${pdma_src_dir}/*.cxx" "${pdma_src_dir}/*.cpp" "${pdma_src_dir}/*.c"
  )
  list(APPEND all_src_files "${pdma_sources}")

  file(GLOB_RECURSE pdmaf_sources
    CONFIGURE_DEPENDS "${pdma_src_dir}/*.F90"
  )
  list(APPEND all_src_fortran_files "${pdmaf_sources}")

  get_all_subdirectories(PDMA_SUBDIRS ${pdma_src_dir})
endif()

set(src_files ${all_src_files})
list(FILTER src_files EXCLUDE REGEX ".*\\.test\\.cpp$")

set(LIB_TYPES "")
if (PDM_ENABLE_STATIC)
  list(APPEND LIB_TYPES STATIC)
endif()
if (PDM_ENABLE_SHARED)
  list(APPEND LIB_TYPES SHARED)
endif()

# if (PDM_ENABLE_CUDA)
#   add_library(pdm_gpu STATIC ${SOURCESGPU})
#   set_target_properties(pdm_gpu PROPERTIES
#                             CUDA_SEPARABLE_COMPILATION ON)
#   target_compile_options(pdm_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#   -maxrregcount=32;--fmad=false;#-G;-g;
#   -gencode=arch=compute_61,code=sm_61;
#   >)
# endif()

add_library(pdm OBJECT ${src_files})

target_include_directories(pdm PUBLIC
  $<BUILD_INTERFACE:${src_dir}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)

foreach(subdir IN LISTS PDM_SUBDIRS)
  target_include_directories(pdm PUBLIC
    $<BUILD_INTERFACE:${subdir}>
  )
endforeach()

foreach(subdir IN LISTS PDMA_SUBDIRS)
  target_include_directories(pdm PUBLIC
    $<BUILD_INTERFACE:${subdir}>
  )
endforeach()

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(pdm PRIVATE ${MATH_LIBRARY})
endif()

target_link_libraries(pdm
PUBLIC
  MPI::MPI_CXX
)

if (PTSCOTCH_FOUND)
  target_link_libraries(pdm
  PUBLIC
    PTScotch::PTScotch
    Scotch::Scotch
    Scotch::ScotchErr
  )
endif()

if (PARMETIS_FOUND)
  target_link_libraries(pdm
  PUBLIC
    ParMETIS::ParMETIS
  )
endif()

if (BLAS_FOUND)
  target_link_libraries(pdm
  PUBLIC
    BLAS::BLAS
  )
endif()

if (LAPACK_FOUND)
  target_link_libraries(pdm
  PUBLIC
    LAPACK::LAPACK
  )
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
get_target_property(PDM_LIB pdm LINK_LIBRARIES)

if (PDM_ENABLE_Fortran)
  if (NOT TARGET MPI::MPI_Fortran)
    find_package(MPI REQUIRED COMPONENTS Fortran)
  endif()
  add_library(pdmf OBJECT ${all_src_fortran_files})
  set_target_properties(pdmf PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_include_directories(pdmf PUBLIC ${PDM_INC})
  target_include_directories(pdmf PUBLIC
    $<INSTALL_INTERFACE:include/modules>
  )
  target_link_libraries(pdmf
  PUBLIC
    MPI::MPI_Fortran
  )
  get_target_property(PDMF_INC pdmf INCLUDE_DIRECTORIES)
  get_target_property(PDMF_LIB pdmf LINK_LIBRARIES)
endif()

foreach(lib_type IN LISTS LIB_TYPES)

  string(TOLOWER "${lib_type}" lib_type_lower)

  add_library(pdm_${lib_type_lower} ${lib_type} $<TARGET_OBJECTS:pdm>)
  add_dependencies(pdm_${lib_type_lower} pdm)
  set_target_properties(pdm_${lib_type_lower} PROPERTIES OUTPUT_NAME "pdm")

  target_include_directories(pdm_${lib_type_lower} PUBLIC ${PDM_INC})
  target_link_libraries(pdm_${lib_type_lower}      PUBLIC ${PDM_LIB})

  if (PDM_ENABLE_Fortran)
    add_library(pdmf_${lib_type_lower} ${lib_type} $<TARGET_OBJECTS:pdmf>)
    add_dependencies(pdmf_${lib_type_lower} pdmf)
    set_target_properties(pdmf_${lib_type_lower} PROPERTIES OUTPUT_NAME "pdmf")

    target_include_directories(pdmf_${lib_type_lower} PUBLIC ${PDMF_INC})
    target_link_libraries(pdmf_${lib_type_lower}      PUBLIC ${PDMF_LIB})

    add_dependencies(pdmf_${lib_type_lower} pdm_${lib_type_lower})
    target_link_libraries(pdmf_${lib_type_lower}
      PUBLIC
        pdm_${lib_type_lower}
    )
  endif()

endforeach()

if (PDM_ENABLE_PYTHON_BINDINGS)
  add_subdirectory(Cython)
endif()

# ----------------------------------------------------------------------
# Unit test
# ----------------------------------------------------------------------
if (PDM_ENABLE_UNIT_TEST)
  include(FetchContent)
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
  )
  FetchContent_GetProperties(doctest)
  if (NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_subdirectory(${doctest_SOURCE_DIR} ${doctest_BINARY_DIR})
  endif()
endif()

# ----------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------
include(CTest)
if (PDM_ENABLE_TESTS)
  add_subdirectory(test)
  if (PDM_ENABLE_EXTENSION_PDMA)
    add_subdirectory(extensions/paradigma/test)
  endif()
endif()

# Configure pdm_run
string(REPLACE ";" " " PDM_LIST_TEST  "${PDM_LIST_TEST}")
string(REPLACE ";" " " PDM_LIST_NRANK "${PDM_LIST_NRANK}")

set (PDM_RUN_PRELOAD "")
set (PDM_RUN_LSAN "")
set (PDM_RUN_UBSAN "")
set (PDM_RUN_PYTHONPATH "")

if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PRELOAD_ASAN OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PDM_RUN_PRELOAD "LD_PRELOAD ${PRELOAD_ASAN}:${PROJECT_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
  set(PDM_RUN_PRELOAD_SRUN "LD_PRELOAD=${PRELOAD_ASAN}:${PROJECT_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
  set(PDM_RUN_LSAN "LSAN_OPTIONS=suppressions=${PROJECT_SOURCE_DIR}/script/asan/asan.supp")
  set(PDM_RUN_UBSAN "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")

  if (PDM_ENABLE_DOC)
    message(FATAL_ERROR "Documention can not be built in Sanitize CMAKE_BUILD_TYPE")
  endif()
endif()

if (DEFINED ENV{PYTHONPATH})
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
else()
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython")
endif()

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_run.in "${PROJECT_BINARY_DIR}/pdm_run"
               @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_launch.in "${PROJECT_BINARY_DIR}/pdm_launch"
               @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_launch_python.in "${PROJECT_BINARY_DIR}/pdm_launch_python"
               @ONLY)

if (PDM_ENABLE_TRAINING)
  add_subdirectory(training)
endif()

# ----------------------------------------------------------------------
# Documentation
# ----------------------------------------------------------------------
add_subdirectory(doc)

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()

# ----------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------
set (_targets "")
if (PDM_ENABLE_Fortran)
  if (PDM_ENABLE_STATIC)
    list (APPEND _targets pdmf_static)
  endif()
  if (PDM_ENABLE_SHARED)
    list (APPEND _targets pdmf_shared)
  endif()
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdm_static)
endif()

if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdm_shared)
endif()

foreach(_target IN ITEMS ${_targets})
  install(TARGETS ${_target} EXPORT ${_target}Targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    ARCHIVE DESTINATION lib
  )
  install(EXPORT ${_target}Targets
    FILE ${_target}Targets.cmake
    NAMESPACE pdm::
    DESTINATION lib/cmake/${_target}
  )
  add_library(pdm::${_target} ALIAS ${_target})
endforeach()

if(PDM_ENABLE_EXTENSION_PDMA)
  file(GLOB_RECURSE all_headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/extensions/paradigma/src/*.h" )
else()
  file(GLOB_RECURSE all_headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
endif()
set(install_headers ${all_headers})
list(FILTER install_headers EXCLUDE REGEX ".*\\_priv\\.h$")

if (PDM_ENABLE_Fortran)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules DESTINATION include)
endif()

install(FILES ${install_headers} DESTINATION include)
install(FILES ${install_configs} DESTINATION include)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/pdmConfig.cmake.in"
  pdmConfig.cmake
  @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdmConfig.cmake"
  DESTINATION lib/cmake/pdm
)


install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindPTScotch.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindParMETIS.cmake
  DESTINATION lib/cmake )


configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/pdmConfigVersion.cmake.in"
  pdmConfigVersion.cmake
  @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdmConfigVersion.cmake"
  DESTINATION lib/cmake/pdm
)

include(pdm_generate_build_defs)
