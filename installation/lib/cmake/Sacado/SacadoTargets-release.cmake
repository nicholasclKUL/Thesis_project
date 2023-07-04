#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Sacado::sacado" for configuration "RELEASE"
set_property(TARGET Sacado::sacado APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Sacado::sacado PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsacado.a"
  )

list(APPEND _cmake_import_check_targets Sacado::sacado )
list(APPEND _cmake_import_check_files_for_Sacado::sacado "${_IMPORT_PREFIX}/lib/libsacado.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
