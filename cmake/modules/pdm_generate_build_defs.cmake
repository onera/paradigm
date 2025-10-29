file(READ ${CMAKE_CURRENT_SOURCE_DIR}/pdm_Build.defs.in pdm_build_defs_file)

# TODO modifier fichier en fonction des options cmake
if (PDM_HAVE_PARMETIS)
  string(REPLACE "#PARMETIS_LIBRARIES" "PARMETIS_LIBRARIES            = ${PARMETIS_LIBRARY}\nMETIS_LIBRARIES               = ${METIS_LIBRARY}"
         pdm_build_defs_file ${pdm_build_defs_file})
else()
  string(REPLACE "#PARMETIS_LIBRARIES" "" pdm_build_defs_file ${pdm_build_defs_file})
endif()

if (PDM_HAVE_PTSCOTCH)
  string(REPLACE "#PTSCOTCH_LIBRARIES" "PTSCOTCH_LIBRARIES            = ${PTScotch_LIBRARY}" pdm_build_defs_file ${pdm_build_defs_file})
else()
  string(REPLACE "#PTSCOTCH_LIBRARIES" "" pdm_build_defs_file ${pdm_build_defs_file})
endif()

if (MPI_C_COMPILER OR MPI_CXX_COMPILER OR MPI_Fortran_COMPILER)
  set(wrapper_mpi "# MPI Wrapper used\n#-----------------------------------------------------------------------")
  if (MPI_C_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_C       = ${MPI_C_COMPILER}")
  endif()
  if (MPI_CXX_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_CXX     = ${MPI_CXX_COMPILER}")
  endif()
  if (MPI_Fortran_COMPILER)
    set(wrapper_mpi "${wrapper_mpi}\nMPI_Fortran = ${MPI_Fortran_COMPILER}")
  endif()
  string(REPLACE "#MPI_Wrapper" ${wrapper_mpi} pdm_build_defs_file ${pdm_build_defs_file})
else()
  string(REPLACE "#MPI_Wrapper" "" pdm_build_defs_file ${pdm_build_defs_file})
endif ()

file(WRITE ${CMAKE_BINARY_DIR}/pdm_Build.defs.in "${pdm_build_defs_file}")

configure_file(${CMAKE_BINARY_DIR}/pdm_Build.defs.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_Build.defs")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pdm_Build.defs
        DESTINATION include)
