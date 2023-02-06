cmake_minimum_required(VERSION 3.20)
cmake_policy (SET CMP0091 NEW)
project (Thesis VERSION 0.0.1)

# Compiler warnings
option(WARNINGS_AS_ERRORS "Enable -Werror or /WX"        Off)
option(ENABLE_FANALYZER   "Enable GCC's -fanalyzer flag" Off)
include(cmake/Warnings.cmake)
add_warnings_target(warnings ${WARNINGS_AS_ERRORS} ${ENABLE_FANALYZER})
add_library(thesis::warnings ALIAS warnings)

if (MSVC)
    set(CMAKE_CXX_FLAGS "/utf-8 /Dnot=! /D__PRETTY_FUNCTION__=__FUNCSIG__ ${CMAKE_CXX_FLAGS}")
endif()

find_package(Taskflow REQUIRED)
find_package(Kokkos REQUIRED)

add_subdirectory(src)
add_subdirectory(applications)