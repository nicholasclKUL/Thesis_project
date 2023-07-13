set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# Set the path to the GCC 8 compiler
set(CMAKE_C_COMPILER /usr/local/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++)

# Set the root directory of the GCC 8 installation (optional)
set(CMAKE_FIND_ROOT_PATH /usr/local/bin/gcc)

# Set other necessary flags or variables specific to GCC 8
# set(CMAKE_CXX_FLAGS "-std=c++20")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
