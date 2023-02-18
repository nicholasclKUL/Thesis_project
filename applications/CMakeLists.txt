add_executable(test #example!!
    "test.cpp" #example!!
)

find_package(alpaqa REQUIRED)
target_link_libraries(test PUBLIC thesis PRIVATE thesis::warnings)
target_link_libraries(test PUBLIC alpaqa::alpaqa alpaqa::warnings)