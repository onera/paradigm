# FindParMETIS.cmake
#
# Définit la target moderne ParMETIS::ParMETIS (et Metis::Metis si possible)
# Utilise PARMETIS_ROOT, PARMETIS_DIR ou $ENV{PARMETIS_ROOT} pour la détection.

set(_PARMETIS_ROOT "")
if(PARMETIS_ROOT)
    set(_PARMETIS_ROOT "${PARMETIS_ROOT}")
elseif(PARMETIS_DIR)
    set(_PARMETIS_ROOT "${PARMETIS_DIR}")
elseif(DEFINED ENV{PARMETIS_DIR})
    set(_PARMETIS_ROOT "$ENV{PARMETIS_DIR}")
elseif(DEFINED ENV{PARMETIS_ROOT})
    set(_PARMETIS_ROOT "$ENV{PARMETIS_ROOT}")
endif()

if(_PARMETIS_ROOT)
    set(_PARMETIS_INCLUDE_HINTS "${_PARMETIS_ROOT}/include")
    set(_PARMETIS_LIB_HINTS "${_PARMETIS_ROOT}/lib" "${_PARMETIS_ROOT}/lib64")
else()
    set(_PARMETIS_INCLUDE_HINTS "")
    set(_PARMETIS_LIB_HINTS "")
endif()

# Recherche des includes
find_path(PARMETIS_INCLUDE_DIR
    NAMES parmetis.h
    HINTS ${_PARMETIS_INCLUDE_HINTS}
)

# Recherche de la bibliothèque ParMETIS
find_library(PARMETIS_LIBRARY
    NAMES parmetis
    HINTS ${_PARMETIS_LIB_HINTS}
)

set(_METIS_ROOT "")
if(METIS_ROOT)
    set(_METIS_ROOT "${METIS_ROOT}")
elseif(METIS_DIR)
    set(_METIS_ROOT "${METIS_DIR}")
elseif(DEFINED ENV{METIS_DIR})
    set(_METIS_ROOT "$ENV{METIS_DIR}")
elseif(DEFINED ENV{METIS_ROOT})
    set(_METIS_ROOT "$ENV{METIS_ROOT}")
endif()

if(NOT _METIS_ROOT)
  # Try to fall back on parmetis install (that could have metis.h)
  if(DEFINED ENV{PARMETIS_DIR})
    set(_METIS_ROOT "$ENV{PARMETIS_DIR}")
  elseif(DEFINED ENV{PARMETIS_ROOT})
    set(_METIS_ROOT "$ENV{PARMETIS_ROOT}")
  endif()
endif()

if(_METIS_ROOT)
    set(_METIS_INCLUDE_HINTS "${_METIS_ROOT}/include")
    set(_METIS_LIB_HINTS "${_METIS_ROOT}/lib" "${_METIS_ROOT}/lib64")
else()
    set(_METIS_INCLUDE_HINTS "")
    set(_METIS_LIB_HINTS "")
endif()

# Recherche des includes
find_path(METIS_INCLUDE_DIR
    NAMES metis.h
    HINTS ${_METIS_INCLUDE_HINTS}
)

# Recherche de la bibliothèque METIS
find_library(METIS_LIBRARY
    NAMES metis
    HINTS ${_METIS_LIB_HINTS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS
    REQUIRED_VARS PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY METIS_LIBRARY
)

if(PARMETIS_FOUND)
    set(PARMETIS_INCLUDE_DIRS "${PARMETIS_INCLUDE_DIR}")
    set(PARMETIS_LIBRARIES "${PARMETIS_LIBRARY};${METIS_LIBRARY}")

    if(NOT TARGET ParMETIS::ParMETIS)
        add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
        set_target_properties(ParMETIS::ParMETIS PROPERTIES
            IMPORTED_LOCATION "${PARMETIS_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "Metis::Metis"
        )
    endif()

    if(NOT TARGET Metis::Metis)
        add_library(Metis::Metis UNKNOWN IMPORTED)
        set_target_properties(Metis::Metis PROPERTIES
            IMPORTED_LOCATION "${METIS_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR}"
        )
    endif()
endif()
