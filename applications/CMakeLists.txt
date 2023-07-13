add_executable(tests #example!!
    "tests.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)

target_link_libraries(tests PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(tests PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(tests PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(tests PUBLIC Kokkos::kokkos)
target_link_libraries(tests PUBLIC Sacado::sacado)