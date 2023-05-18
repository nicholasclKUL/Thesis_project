add_executable(test_alm_AD #example!!
    "test_alm_AD.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)

target_link_libraries(test_alm_AD PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(test_alm_AD PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(test_alm_AD PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(test_alm_AD PUBLIC Kokkos::kokkos)
target_link_libraries(test_alm_AD PUBLIC Sacado::sacado)
