# --------------------------
# Construction des Tests
# --------------------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif() 

if (MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif() 

unset(LINK_LIBRARIES)

if (PDM_ENABLE_STATIC)
 list(APPEND LINK_LIBRARIES
      pdm_static
      pdm_mpi_static
      ${COMMON_LIBRARIES})
else ()
 list(APPEND LINK_LIBRARIES
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)                          


# Test part_cube
# --------------

file(GLOB_RECURSE SOURCES 
     part_cube.c)

add_executable(part_cube ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(part_cube
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(part_cube
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(part_cube PRIVATE ${CMAKE_SOURCE_DIR}
                                     PRIVATE ${CMAKE_BINARY_DIR}
                                     PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(part_cube PRIVATE ${PDM_INC})

target_link_libraries(part_cube ${LINK_LIBRARIES})

install(TARGETS part_cube
        RUNTIME DESTINATION bin)

add_test (part_cube ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/part_cube
   ${MPIEXEC_POSTFLAGS})

message(STATUS "Add test part_cube - Succes")

# Test 3D (plusieurs sous-domaines par proc)
# ------------------------------------------

unset(SOURCES)
file(GLOB_RECURSE SOURCES 
     cs_test_3d.c)

add_executable(cs_test_3d ${SOURCES})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(cs_test_3d
    PROPERTIES
    COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_C_COMPILE_FLAGS)
  set_target_properties(cs_test_3d
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

target_include_directories(cs_test_3d PRIVATE ${CMAKE_SOURCE_DIR}
                                      PRIVATE ${CMAKE_BINARY_DIR}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR/../src/msg}
                                      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(cs_test_3d PRIVATE ${PDM_INC})
target_link_libraries(cs_test_3d ${LINK_LIBRARIES})

install(TARGETS cs_test_3d
        RUNTIME DESTINATION bin)

# Command to launch by CTest.

add_test (cs_test_3d ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
   ${MPIEXEC_PREFLAGS}
   ${CMAKE_CURRENT_BINARY_DIR}/cs_test_3d
   ${MPIEXEC_POSTFLAGS})
