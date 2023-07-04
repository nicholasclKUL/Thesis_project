add_executable(tests_ss_vs_ms #example!!
    "tests_ss_vs_ms.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)

target_link_libraries(tests_ss_vs_ms PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(tests_ss_vs_ms PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(tests_ss_vs_ms PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(tests_ss_vs_ms PUBLIC Kokkos::kokkos)
target_link_libraries(tests_ss_vs_ms PUBLIC Sacado::sacado)
