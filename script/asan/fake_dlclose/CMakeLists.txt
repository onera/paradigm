
add_library(dlclose SHARED dlclose.c)

target_include_directories(dlclose
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../../src>
)
 
