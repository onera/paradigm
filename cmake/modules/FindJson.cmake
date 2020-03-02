# Find the json library
#
# JSON_FOUND - True if json was found.
# JSON_LIBRARIES - The libraries needed to use json
# JSON_INCLUDE_DIRS - Location of json.h

set(JSON_FOUND FALSE)
set(JSON_ROOT  /stck/jcoulet/dev/dev-Test/json-c/install)

find_path(JSON_INCLUDE_DIRS json.h
	HINTS ${JSON_ROOT}
	PATH_SUFFIXES include/json-c
  DOC "Directory where the json header is located"
  )

find_library(JSON_LIBRARIES json-c
	HINTS ${JSON_ROOT}
	PATH_SUFFIXES lib64
  NO_DEFAULT_PATH
  DOC "The json header library"
  )

if (JSON_LIBRARIES)
	string(REGEX REPLACE "(^.*)/json-c.*$" "\\1" JSON_LIBRARY_PATH ${JSON_LIBRARIES} )
endif (JSON_LIBRARIES)

IF(JSON_FIND_REQUIRED AND NOT JSON_FOUND)
  message(SEND_ERROR "Unable to find the requested Json libraries.")
ENDIF(JSON_FIND_REQUIRED AND NOT JSON_FOUND)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(JSON DEFAULT_MSG JSON_LIBRARIES JSON_INCLUDE_DIRS JSON_LIBRARY_PATH)

MARK_AS_ADVANCED(
	JSON_INCLUDE_DIRS
  JSON_LIBRARIES
  JSON_LIBRARY_PATH
)

