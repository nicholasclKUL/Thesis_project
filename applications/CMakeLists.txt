add_executable(test2 #example!!
    "test2.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)

target_link_libraries(test2 PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(test2 PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(test2 PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(test2 PUBLIC Kokkos::kokkos)
