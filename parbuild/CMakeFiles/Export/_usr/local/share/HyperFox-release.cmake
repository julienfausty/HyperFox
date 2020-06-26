#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hyperfox" for configuration "Release"
set_property(TARGET hyperfox APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hyperfox PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhyperfox.so"
  IMPORTED_SONAME_RELEASE "libhyperfox.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS hyperfox )
list(APPEND _IMPORT_CHECK_FILES_FOR_hyperfox "${_IMPORT_PREFIX}/lib/libhyperfox.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
