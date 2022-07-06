#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libtommath" for configuration "Release"
set_property(TARGET libtommath APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libtommath PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtommath.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS libtommath )
list(APPEND _IMPORT_CHECK_FILES_FOR_libtommath "${_IMPORT_PREFIX}/lib/libtommath.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
