add_executable(test_smallhangingchaindynamics #example!!
    "test_smallhangingchaindynamics.cpp" #example!!
)

find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)

target_link_libraries(test_smallhangingchaindynamics PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(test_smallhangingchaindynamics PUBLIC problems PRIVATE problems::warnings)
target_link_libraries(test_smallhangingchaindynamics PUBLIC alpaqa::alpaqa alpaqa::warnings)
target_link_libraries(test_smallhangingchaindynamics PUBLIC Kokkos::kokkos)
target_link_libraries(test_smallhangingchaindynamics PUBLIC Sacado::sacado)
