include(UseCython)

if (PDM_ENABLE_EXTENSION_PDMA)
  file(GLOB SOURCE_PXI
       CONFIGURE_DEPENDS 
       *.pxi
       ${PDM_SOURCE_DIR}/extensions/paradigma/Cython/*.pxi)

  file(GLOB SOURCE_PY
       CONFIGURE_DEPENDS 
       *.py
       ${PDM_SOURCE_DIR}/extensions/paradigma/Cython/*.py)
else ()
  file(GLOB SOURCE_PXI
       CONFIGURE_DEPENDS 
       *.pxi)

  file(GLOB SOURCE_PY
       CONFIGURE_DEPENDS 
       *.py)
endif()

foreach(F IN LISTS SOURCE_PXI SOURCE_PY)
    get_filename_component(F_BASENAME "${F}" NAME )    
    configure_file(${F} ${CMAKE_CURRENT_BINARY_DIR}/${F_BASENAME} ) 
endforeach()

configure_file(Pypdm.pyx.in Pypdm.pyx)
configure_file(pdm_part_coarse.pxi.in pdm_part_coarse.pxi)

set(_pyx_file "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.pyx")
set_source_files_properties(${_pyx_file} PROPERTIES CYTHON_IS_CXX FALSE)
# set_source_files_properties(${_pyx_file} PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module (Pypdm
                   PYX_SOURCES "${_pyx_file}"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI} pdm_SHARED)

set_target_properties(Pypdm PROPERTIES PREFIX "") # do not prefix by lib

set_target_properties(Pypdm PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/")

target_link_libraries(Pypdm PRIVATE pdm_SHARED)

set (SITE_PACKAGES_PDM "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/Pypdm")

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "Clang")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (Pypdm PROPERTIES COMPILE_FLAGS "-w")
endif()

install(DIRECTORY DESTINATION ${SITE_PACKAGES_PDM})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_PDM}")

install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.so"
        DESTINATION            "${SITE_PACKAGES_PDM}")
