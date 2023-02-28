add_library(thesis
"src/Thesis.cpp"

"include/Thesis/Thesis.hpp"
)

target_link_libraries(thesis PRIVATE Kokkos::kokkos)
target_link_libraries(thesis PRIVATE Eigen3::Eigen)
target_link_libraries(thesis PRIVATE alpaqa::alpaqa)

target_link_libraries(thesis PRIVATE thesis::warnings)
target_include_directories(thesis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
