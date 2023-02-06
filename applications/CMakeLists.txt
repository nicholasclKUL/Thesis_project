add_executable(hello #example!!
    "hello.cpp" #example!!
)
target_link_libraries(hello PUBLIC thesis PRIVATE thesis::warnings)