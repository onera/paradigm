file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/prepare_ci DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  add_subdirectory(asan)	
endif()

set (PDM_RUN_PRELOAD "")
set (PDM_RUN_LSAN "")
set (PDM_RUN_PYTHONPATH "")

if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PRELOAD_ASAN OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PDM_RUN_PRELOAD "LD_PRELOAD=\"${PRELOAD_ASAN} ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so\"")

  set(PDM_RUN_LSAN "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/script/asan/asan.supp")
endif()

if(DEFINED ENV{PYTHONPATH})
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
else()
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${CMAKE_BINARY_DIR}/Cython")
endif()

configure_file(pdm_run.in "${CMAKE_BINARY_DIR}/test/pdm_run")
file(CHMOD ${CMAKE_BINARY_DIR}/test/pdm_run
  PERMISSIONS  OWNER_READ OWNER_EXECUTE)
if (PDM_ENABLE_EXTENSION_PDMA)
  configure_file(pdm_run.in "${CMAKE_BINARY_DIR}/extension/paradigma/test/pdm_run")
  file(CHMOD ${CMAKE_BINARY_DIR}/extension/paradigma/test/pdm_run
       PERMISSIONS  OWNER_READ OWNER_EXECUTE)
endif()