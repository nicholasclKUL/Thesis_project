add_library(thesis
    "src/Thesis.cpp" #example!!
    "include/Thesis/Thesis.hpp" #example!!
)

target_link_libraries(thesis PRIVATE Taskflow::Taskflow)
target_link_libraries(thesis PRIVATE Kokkos::kokkos)

target_link_libraries(thesis PRIVATE thesis::warnings)
target_include_directories(thesis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)