add_executable(OCP #example!!
    "OCP.cpp" #example!!
)

find_package(Eigen3 REQUIRED)
find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)

target_link_libraries(OCP PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(OCP PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(OCP PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(OCP PUBLIC Kokkos::kokkos)
target_link_libraries(OCP PUBLIC Sacado::sacado)