################################################################################
#
# test_c_create and test_fortran_create add a Fortran test or a C test
#
# They uses LINK_LIBRARIES and TEST_INC variables
#
################################################################################


function (add_test_pdm_run name n_proc LIST_TEST LIST_NRANK)

    set (${LIST_TEST} ${${LIST_TEST}} "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    set (${LIST_NRANK} ${${LIST_NRANK}} "${n_proc}")

    set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
    set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

endfunction()

function(test_create names n_procs LIST_TEST LIST_NRANK)

  set (command_test "")

  set (MPIEXEC_SPLIT "")

  set (str_name "")
  set (str_n_proc "")
  set (str_sep "")

  foreach(_name _n_proc IN ZIP_LISTS names n_procs)
    get_filename_component(_name_base ${_name} NAME_WE)
    get_filename_component(_name_ext ${_name} LAST_EXT)

    set (PRE_EXE "")
    set (POST_EXE "")

    set (str_n_proc "${str_n_proc}${str_sep}${_n_proc}")

    if ("${_name_ext}" STREQUAL ".c")
      set (str_name "${str_name}${str_sep}${_name_base}")
      add_executable(${_name_base} ${_name})

      if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
        set_target_properties(${_name_base}
                              PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
      endif()
      target_include_directories(${_name_base} PRIVATE ${CMAKE_SOURCE_DIR}
                                         PRIVATE ${CMAKE_BINARY_DIR}
                                         PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
      target_include_directories(${_name_base} PRIVATE ${TEST_INC})
      target_link_libraries(${_name_base} ${LINK_LIBRARIES})
      if(mmg_FOUND)
       target_link_libraries(${_name_base} Mmg::libmmg2d_so)
      endif()
      if (LAPACK_FOUND)
        target_link_libraries(${_name_base} LAPACK::LAPACK)
      endif()
      if (NOT LAPACK_FOUND AND BLAS_FOUND)
        target_link_libraries(${_name_base} BLAS::BLAS)
      endif()

    elseif (("${_name_ext}" STREQUAL ".f90") OR ("${_name_ext}" STREQUAL ".F90"))
      set (str_name "${str_name}${str_sep}${_name_base}")
      add_executable(${_name_base} ${_name})
      if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
        set_target_properties(${_name_base}
                              PROPERTIES
                              COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
       target_include_directories(${_name_base} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
      endif()
      target_include_directories(${_name_base} PRIVATE ${CMAKE_SOURCE_DIR}
                                         PRIVATE ${CMAKE_BINARY_DIR}
                                         PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
      target_include_directories(${_name_base} PRIVATE ${TEST_INC})
      target_link_libraries(${_name_base} ${LINK_LIBRARIES})
      set_target_properties(${_name_base} PROPERTIES LINKER_LANGUAGE "Fortran")

    elseif ("${_name_ext}" STREQUAL ".py")
      set (str_name "${str_name}${str_sep}${_name_base}")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_name} ${CMAKE_CURRENT_BINARY_DIR}/${_name})

      add_custom_target(${_name}
                          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_name})
      set (PRE_EXE "${Python_EXECUTABLE}")
      set (POST_EXE ".py")
    else ()
      message (FATAL_ERROR "unknown extension file ${_name_ext} ${_name}")
    endif()

    if (NOT test_name)
      set (test_name ${_name_base})
    endif()

    list (APPEND target_exe ${_name_base})
    list (APPEND command_test ${MPIEXEC_SPLIT} ${MPIEXEC_NUMPROC_FLAG} ${_n_proc} ${MPIEXEC_PREFLAGS} ${MPIEXEC_GENV_COMMAND} ${MPIEXEC_GENV_PRELOAD} ${MPIEXEC_GENV_PRELOAD_PATH} ${PRE_EXE} ${CMAKE_CURRENT_BINARY_DIR}/${_name_base}${POST_EXE})
    set (MPIEXEC_SPLIT ":")
    set (str_sep ":")
  endforeach()

  add_test (${test_name} ${MPIEXEC} ${command_test} ${MPIEXEC_POSTFLAGS})

  set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
  set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if(DEFINED ENV{PATH})
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
  else()
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${test_name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()

  add_test_pdm_run (${str_name} ${str_n_proc} ${LIST_TEST} ${LIST_NRANK})

  set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
  set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

endfunction()

function(test_c_create name n_proc LIST_TEST LIST_NRANK)
  add_executable(${name} "${name}.c")

  if(PDM_ENABLE_STATIC AND NOT PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdm_STATIC)
  endif()
  if(PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdm_SHARED)
  endif()

  install(TARGETS ${name} RUNTIME DESTINATION bin)

  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${MPIEXEC_GENV_COMMAND}
            ${MPIEXEC_GENV_PRELOAD}
            ${MPIEXEC_GENV_PRELOAD_PATH}
            ${CMAKE_CURRENT_BINARY_DIR}/${name}
            ${MPIEXEC_POSTFLAGS})

  add_test_pdm_run (${name} ${n_proc} ${LIST_TEST} ${LIST_NRANK})

  set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
  set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    # execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PYTHON_TEST_ENV1 OUTPUT_STRIP_TRAILING_WHITESPACE)
    # set(PYTHON_TEST_ENV "LD_PRELOAD=${PYTHON_TEST_ENV1} ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
    # list(APPEND LIST_TEST_ENV "${PYTHON_TEST_ENV}")
    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if(DEFINED ENV{PATH})
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
  else()
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()

endfunction()

function(test_fortran_create name n_proc LIST_TEST LIST_NRANK)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.f90")
    add_executable(${name} "${name}.f90")
  else ()
    add_executable(${name} "${name}.F90")
  endif()

  if(PDM_ENABLE_STATIC AND NOT PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdmf_STATIC)
    get_target_property(MODDIR pdmf Fortran_MODULE_DIRECTORY)
  endif()
  if(PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdmf_SHARED)
    get_target_property(MODDIR pdmf Fortran_MODULE_DIRECTORY)
  endif()

  target_include_directories(${name} PRIVATE ${MODDIR})
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE "Fortran")

  install(TARGETS ${name} RUNTIME DESTINATION bin)
  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${MPIEXEC_GENV_COMMAND}
            ${MPIEXEC_GENV_PRELOAD}
            ${MPIEXEC_GENV_PRELOAD_PATH}
            ${CMAKE_CURRENT_BINARY_DIR}/${name}
            ${MPIEXEC_POSTFLAGS})

   add_test_pdm_run (${name} ${n_proc} ${LIST_TEST} ${LIST_NRANK})
    set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
    set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    # execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PYTHON_TEST_ENV1 OUTPUT_STRIP_TRAILING_WHITESPACE)
    # set(PYTHON_TEST_ENV "LD_PRELOAD=${PYTHON_TEST_ENV1} ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
    # list(APPEND LIST_TEST_ENV "${PYTHON_TEST_ENV}")

    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if(DEFINED ENV{PATH})
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
  else()
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()

endfunction()

function(test_python_create name n_proc LIST_TEST LIST_NRANK)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)

  add_custom_target(${name}
                      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py)

  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${MPIEXEC_GENV_COMMAND}
            ${MPIEXEC_GENV_PRELOAD}
            ${MPIEXEC_GENV_PRELOAD_PATH}
            ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${name}.py
            ${MPIEXEC_POSTFLAGS})

  if (TARGET Pypdm)
    add_dependencies(${name} Pypdm)
  endif ()

  add_test_pdm_run (${name} ${n_proc} ${LIST_TEST} ${LIST_NRANK})
    set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
    set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    # execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PYTHON_TEST_ENV1 OUTPUT_STRIP_TRAILING_WHITESPACE)
    # set(PYTHON_TEST_ENV "LD_PRELOAD=${PYTHON_TEST_ENV1} ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
    # list(APPEND LIST_TEST_ENV "${PYTHON_TEST_ENV}")
    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if(DEFINED ENV{PYTHONPATH})
    list(APPEND LIST_TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
  else()
    list(APPEND LIST_TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython")
  endif()

  if(DEFINED ENV{PATH})
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
  else()
    list(APPEND LIST_TEST_ENV "PATH=${CMAKE_BINARY_DIR}/bin")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()

endfunction()

function(test_cpp_unit_create name n_proc LIST_TEST LIST_NRANK)
  set(options)
  set(one_value_args)
  set(multi_value_args SOURCES)
  cmake_parse_arguments(ARGS "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  add_executable(${name} "${name}.cpp" ${ARGS_SOURCES})

  if(PDM_ENABLE_STATIC AND NOT PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdm_STATIC)
  endif()
  if(PDM_ENABLE_SHARED)
    target_link_libraries(${name} PRIVATE pdm_SHARED)
  endif()

  target_link_libraries(${name} PRIVATE doctest::doctest)
  install(TARGETS ${name} RUNTIME DESTINATION bin)
  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${CMAKE_CURRENT_BINARY_DIR}/${name}
            ${MPIEXEC_POSTFLAGS})

  add_test_pdm_run (${name} ${n_proc} ${LIST_TEST} ${LIST_NRANK})
  set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
  set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)


  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    # execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PYTHON_TEST_ENV1 OUTPUT_STRIP_TRAILING_WHITESPACE)
    # set(PYTHON_TEST_ENV "LD_PRELOAD=${PYTHON_TEST_ENV1} ${PDM_BINARY_DLCLOSE_DIR}/libdlclose.so")
    # list(APPEND LIST_TEST_ENV "${PYTHON_TEST_ENV}")

    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()
endfunction()

function(test_cpp_create name n_proc LIST_TEST LIST_NRANK)
   add_executable(${name} "${name}.cpp")
   if ((NOT MPI_CXX_COMPILER) AND MPI_CXX_COMPILE_FLAGS)
     set_target_properties(${name}
                           PROPERTIES
                           COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
   endif()
   target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
   target_include_directories(${name} PRIVATE ${TEST_INC})
   target_link_libraries(${name} ${LINK_LIBRARIES})
   install(TARGETS ${name} RUNTIME DESTINATION bin)
   add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
             ${MPIEXEC_PREFLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${name}
             ${MPIEXEC_POSTFLAGS})

   add_test_pdm_run (${name} ${n_proc} ${LIST_TEST} ${LIST_NRANK})
   set (${LIST_TEST} ${${LIST_TEST}} PARENT_SCOPE)
   set (${LIST_NRANK} ${${LIST_NRANK}} PARENT_SCOPE)

  set (LIST_TEST_ENV "")

  if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    # execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PYTHON_TEST_ENV1 OUTPUT_STRIP_TRAILING_WHITESPACE)
    # set(PYTHON_TEST_ENV "LD_PRELOAD=${PYTHON_TEST_ENV1} ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
    # list(APPEND LIST_TEST_ENV "${PYTHON_TEST_ENV}")
    list(APPEND LIST_TEST_ENV "LSAN_OPTIONS=suppressions=${PDM_SOURCE_DIR}/script/asan/asan.supp")
    list(APPEND LIST_TEST_ENV "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0")
  endif()

  if (LIST_TEST_ENV)
    set_property(TEST ${name} PROPERTY ENVIRONMENT "${LIST_TEST_ENV}")
  endif()

endfunction()
