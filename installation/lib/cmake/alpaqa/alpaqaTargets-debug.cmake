#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "alpaqa::alpaqa" for configuration "Debug"
set_property(TARGET alpaqa::alpaqa APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::alpaqa PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libalpaqa_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::alpaqa )
list(APPEND _cmake_import_check_files_for_alpaqa::alpaqa "${_IMPORT_PREFIX}/lib/libalpaqa_d.a" )

# Import target "alpaqa::casadi-loader" for configuration "Debug"
set_property(TARGET alpaqa::casadi-loader APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::casadi-loader PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libalpaqa-casadi-loader_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::casadi-loader )
list(APPEND _cmake_import_check_files_for_alpaqa::casadi-loader "${_IMPORT_PREFIX}/lib/libalpaqa-casadi-loader_d.a" )

# Import target "alpaqa::casadi-ocp-loader" for configuration "Debug"
set_property(TARGET alpaqa::casadi-ocp-loader APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::casadi-ocp-loader PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libalpaqa-casadi-ocp-loader_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::casadi-ocp-loader )
list(APPEND _cmake_import_check_files_for_alpaqa::casadi-ocp-loader "${_IMPORT_PREFIX}/lib/libalpaqa-casadi-ocp-loader_d.a" )

# Import target "alpaqa::dl-loader" for configuration "Debug"
set_property(TARGET alpaqa::dl-loader APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::dl-loader PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libalpaqa-dl-loader_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::dl-loader )
list(APPEND _cmake_import_check_files_for_alpaqa::dl-loader "${_IMPORT_PREFIX}/lib/libalpaqa-dl-loader_d.a" )

# Import target "alpaqa::lbfgsb-fortran" for configuration "Debug"
set_property(TARGET alpaqa::lbfgsb-fortran APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::lbfgsb-fortran PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "Fortran"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/liblbfgsb-fortran_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::lbfgsb-fortran )
list(APPEND _cmake_import_check_files_for_alpaqa::lbfgsb-fortran "${_IMPORT_PREFIX}/lib/liblbfgsb-fortran_d.a" )

# Import target "alpaqa::lbfgsb-adapter" for configuration "Debug"
set_property(TARGET alpaqa::lbfgsb-adapter APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::lbfgsb-adapter PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX;Fortran"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libalpaqa-lbfgsb-adapter_d.a"
  )

list(APPEND _cmake_import_check_targets alpaqa::lbfgsb-adapter )
list(APPEND _cmake_import_check_files_for_alpaqa::lbfgsb-adapter "${_IMPORT_PREFIX}/lib/libalpaqa-lbfgsb-adapter_d.a" )

# Import target "alpaqa::driver" for configuration "Debug"
set_property(TARGET alpaqa::driver APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(alpaqa::driver PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/alpaqa-driver_d"
  )

list(APPEND _cmake_import_check_targets alpaqa::driver )
list(APPEND _cmake_import_check_files_for_alpaqa::driver "${_IMPORT_PREFIX}/bin/alpaqa-driver_d" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
