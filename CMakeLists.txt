cmake_minimum_required(VERSION 3.20)
cmake_policy (SET CMP0091 NEW)

# Set toolchain for cross-compilation
# set(CMAKE_TOOLCHAIN_FILE /home/nicholas/toolchains/x86_64-centos7-linux-gnu/cmake/x86_64-centos7-linux-gnu.toolchain.cmake)
add_compile_options(-finput-charset=UTF-8)

project (Thesis VERSION 0.0.1)

# Compiler warnings
option(WARNINGS_AS_ERRORS "Enable -Werror or /WX"        Off)
option(ENABLE_FANALYZER   "Enable GCC's -fanalyzer flag" Off)
include(cmake/Warnings.cmake)
add_warnings_target(warnings ${WARNINGS_AS_ERRORS} ${ENABLE_FANALYZER})
add_library(thesis::warnings ALIAS warnings)
add_library(problems::warnings ALIAS warnings)

if (MSVC)
    set(CMAKE_CXX_FLAGS "/utf-8 /Dnot=! /D__PRETTY_FUNCTION__=__FUNCSIG__ ${CMAKE_CXX_FLAGS}")
endif()

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

set(CMAKE_POSITION_INDEPENDENT_CODE On)

find_package(OpenMP REQUIRED) #add OpenMP::OpenMP<CXX> to the file I am calling it
find_package(Eigen3 REQUIRED)
find_package(alpaqa REQUIRED)
find_package(Kokkos REQUIRED)
find_package(Sacado REQUIRED)
find_package(Threads REQUIRED)

add_subdirectory(src)
add_subdirectory(problems)
add_subdirectory(applications)