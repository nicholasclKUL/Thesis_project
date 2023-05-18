add_library(problems

# benchmark problems
"linear_dynamics.hpp"
"nonlinear_example1.hpp"
"hanging_chain.hpp"
"quadcopter.hpp"

# benchmark problems using AD
"quadcopter_AD.hpp"
)

target_link_libraries(problems PRIVATE Kokkos::kokkos)
target_link_libraries(problems PRIVATE Eigen3::Eigen)
target_link_libraries(problems PRIVATE alpaqa::alpaqa)
target_link_libraries(problems PUBLIC thesis PRIVATE thesis::warnings)

target_link_libraries(problems PRIVATE problems::warnings)
target_include_directories(problems PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)