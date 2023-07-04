#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "KokkosCore::kokkoscore" for configuration "RELEASE"
set_property(TARGET KokkosCore::kokkoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(KokkosCore::kokkoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkokkoscore.a"
  )

list(APPEND _cmake_import_check_targets KokkosCore::kokkoscore )
list(APPEND _cmake_import_check_files_for_KokkosCore::kokkoscore "${_IMPORT_PREFIX}/lib/libkokkoscore.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
