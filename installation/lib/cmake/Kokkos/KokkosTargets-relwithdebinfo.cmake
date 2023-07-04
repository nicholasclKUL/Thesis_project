#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Kokkos::kokkoscore" for configuration "RelWithDebInfo"
set_property(TARGET Kokkos::kokkoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Kokkos::kokkoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libkokkoscore.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkoscore )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkoscore "${_IMPORT_PREFIX}/lib/libkokkoscore.a" )

# Import target "Kokkos::kokkoscontainers" for configuration "RelWithDebInfo"
set_property(TARGET Kokkos::kokkoscontainers APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Kokkos::kokkoscontainers PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libkokkoscontainers.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkoscontainers )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkoscontainers "${_IMPORT_PREFIX}/lib/libkokkoscontainers.a" )

# Import target "Kokkos::kokkossimd" for configuration "RelWithDebInfo"
set_property(TARGET Kokkos::kokkossimd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Kokkos::kokkossimd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libkokkossimd.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkossimd )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkossimd "${_IMPORT_PREFIX}/lib/libkokkossimd.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
