# FindPTScotch.cmake
#
# Définit les targets modernes pour toutes les librairies PT-Scotch et Scotch

# Gestion du chemin d'installation personnalisé
set(_PTSCOTCH_ROOT "")
if(PTSCOTCH_ROOT)
    set(_PTSCOTCH_ROOT "${PTSCOTCH_ROOT}")
elseif(DEFINED ENV{PTSCOTCH_ROOT})
    set(_PTSCOTCH_ROOT "$ENV{PTSCOTCH_ROOT}")
endif()

if(_PTSCOTCH_ROOT)
    set(_PTSCOTCH_INCLUDE_HINTS "${_PTSCOTCH_ROOT}/include")
    set(_PTSCOTCH_LIB_HINTS "${_PTSCOTCH_ROOT}/lib" "${_PTSCOTCH_ROOT}/lib64")
else()
    set(_PTSCOTCH_INCLUDE_HINTS "")
    set(_PTSCOTCH_LIB_HINTS "")
endif()

# Recherche des includes
find_path(PTScotch_INCLUDE_DIR
    NAMES ptscotch.h
    HINTS ${_PTSCOTCH_INCLUDE_HINTS}
)

# Recherche des librairies
macro(_find_ptscotch_lib _var _name)
    find_library(${_var}
        NAMES ${_name}
        HINTS ${_PTSCOTCH_LIB_HINTS}
    )
endmacro()

_find_ptscotch_lib(PTScotch_LIBRARY ptscotch)
_find_ptscotch_lib(PTScotchErr_LIBRARY ptscotcherr)
_find_ptscotch_lib(PTScotchErrExit_LIBRARY ptscotcherrexit)
_find_ptscotch_lib(PTScotchParmetis_LIBRARY ptscotchparmetis)
_find_ptscotch_lib(PTESMumps_LIBRARY ptesmumps)
_find_ptscotch_lib(ESMumps_LIBRARY esmumps)

_find_ptscotch_lib(Scotch_LIBRARY scotch)
_find_ptscotch_lib(ScotchErr_LIBRARY scotcherr)
_find_ptscotch_lib(ScotchErrExit_LIBRARY scotcherrexit)
_find_ptscotch_lib(ScotchMetis_LIBRARY scotchmetis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTScotch
    REQUIRED_VARS PTScotch_INCLUDE_DIR PTScotch_LIBRARY
)

if(PTScotch_FOUND)
    set(PTScotch_INCLUDE_DIRS "${PTScotch_INCLUDE_DIR}")

    # Targets pour PT-Scotch
    if(NOT TARGET PTScotch::PTScotch)
        add_library(PTScotch::PTScotch UNKNOWN IMPORTED)
        set_target_properties(PTScotch::PTScotch PROPERTIES
            IMPORTED_LOCATION "${PTScotch_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(PTScotchErr_LIBRARY AND NOT TARGET PTScotch::PTScotchErr)
        add_library(PTScotch::PTScotchErr UNKNOWN IMPORTED)
        set_target_properties(PTScotch::PTScotchErr PROPERTIES
            IMPORTED_LOCATION "${PTScotchErr_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(PTScotchErrExit_LIBRARY AND NOT TARGET PTScotch::PTScotchErrExit)
        add_library(PTScotch::PTScotchErrExit UNKNOWN IMPORTED)
        set_target_properties(PTScotch::PTScotchErrExit PROPERTIES
            IMPORTED_LOCATION "${PTScotchErrExit_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(PTScotchParmetis_LIBRARY AND NOT TARGET PTScotch::PTScotchParmetis)
        add_library(PTScotch::PTScotchParmetis UNKNOWN IMPORTED)
        set_target_properties(PTScotch::PTScotchParmetis PROPERTIES
            IMPORTED_LOCATION "${PTScotchParmetis_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(PTESMumps_LIBRARY AND NOT TARGET PTScotch::PTESMumps)
        add_library(PTScotch::PTESMumps UNKNOWN IMPORTED)
        set_target_properties(PTScotch::PTESMumps PROPERTIES
            IMPORTED_LOCATION "${PTESMumps_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(ESMumps_LIBRARY AND NOT TARGET PTScotch::ESMumps)
        add_library(PTScotch::ESMumps UNKNOWN IMPORTED)
        set_target_properties(PTScotch::ESMumps PROPERTIES
            IMPORTED_LOCATION "${ESMumps_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    # Targets pour Scotch (séquentiel)
    if(Scotch_LIBRARY AND NOT TARGET Scotch::Scotch)
        add_library(Scotch::Scotch UNKNOWN IMPORTED)
        set_target_properties(Scotch::Scotch PROPERTIES
            IMPORTED_LOCATION "${Scotch_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(ScotchErr_LIBRARY AND NOT TARGET Scotch::ScotchErr)
        add_library(Scotch::ScotchErr UNKNOWN IMPORTED)
        set_target_properties(Scotch::ScotchErr PROPERTIES
            IMPORTED_LOCATION "${ScotchErr_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(ScotchErrExit_LIBRARY AND NOT TARGET Scotch::ScotchErrExit)
        add_library(Scotch::ScotchErrExit UNKNOWN IMPORTED)
        set_target_properties(Scotch::ScotchErrExit PROPERTIES
            IMPORTED_LOCATION "${ScotchErrExit_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()

    if(ScotchMetis_LIBRARY AND NOT TARGET Scotch::ScotchMetis)
        add_library(Scotch::ScotchMetis UNKNOWN IMPORTED)
        set_target_properties(Scotch::ScotchMetis PROPERTIES
            IMPORTED_LOCATION "${ScotchMetis_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}"
        )
    endif()
endif()

