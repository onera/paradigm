set (PDM_DEF_VERSION "2.3.0")
#set (PDM_DEF_VERSION "2.0.0-a1.untagged")
string(REPLACE "-" ";" PDM_DEF_VERSION_LIST ${PDM_DEF_VERSION})
list(GET PDM_DEF_VERSION_LIST 0 PDM_DEF_VERSION_MAJOR)
