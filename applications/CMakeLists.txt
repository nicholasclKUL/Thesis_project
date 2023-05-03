add_executable(tests_alm #example!!
    "tests_alm.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)

target_link_libraries(tests_alm PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(tests_alm PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(tests_alm PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(tests_alm PUBLIC Kokkos::kokkos)
