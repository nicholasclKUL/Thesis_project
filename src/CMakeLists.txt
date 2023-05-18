add_library(thesis 

# solver
"include/thesis/thesis.hpp"
"include/thesis/para-panoc.hpp"
"include/thesis/para-alm.hpp"
#"include/thesis/para-lbfgs.hpp"
"include/thesis/printing.hpp"
#"include/thesis/ms-accelerators.hpp"
"include/thesis/para-panoc-helpers.hpp"
"include/thesis/auto-diff-sacado-helpers.hpp"
)

target_link_libraries(thesis PRIVATE Kokkos::kokkos)
target_link_libraries(thesis PRIVATE Eigen3::Eigen)
target_link_libraries(thesis PRIVATE alpaqa::alpaqa)
target_link_libraries(thesis PRIVATE Sacado::sacado)

target_link_libraries(thesis PRIVATE thesis::warnings)
target_include_directories(thesis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
