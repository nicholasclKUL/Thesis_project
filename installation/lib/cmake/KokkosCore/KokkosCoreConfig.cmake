# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

##############################################################################
#
# CMake variable for use by Trilinos/KokkosCore clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "KokkosCore requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/KokkosCore build
## ---------------------------------------------------------------------------

set(KokkosCore_CXX_COMPILER "/usr/bin/g++")

set(KokkosCore_C_COMPILER "/usr/bin/gcc")

set(KokkosCore_Fortran_COMPILER "/usr/bin/gfortran")
# Deprecated!
set(KokkosCore_FORTRAN_COMPILER "/usr/bin/gfortran") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/KokkosCore build
## ---------------------------------------------------------------------------

## Give the build type
set(KokkosCore_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(KokkosCore_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(KokkosCore_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(KokkosCore_Fortran_FLAGS [[ -O3]])
# Deprecated
set(KokkosCore_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(KokkosCore_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(KokkosCore_SHARED_LIB_RPATH_COMMAND "")
set(KokkosCore_BUILD_SHARED_LIBS "FALSE")

set(KokkosCore_LINKER /usr/bin/ld)
set(KokkosCore_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(KokkosCore_INSTALL_DIR "/home/nicholas/trilinos")

## List of package libraries
set(KokkosCore_LIBRARIES KokkosCore::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(KokkosCore_MPI_LIBRARIES "")
set(KokkosCore_MPI_LIBRARY_DIRS "")
set(KokkosCore_MPI_INCLUDE_DIRS "")
set(KokkosCore_MPI_EXEC "")
set(KokkosCore_MPI_EXEC_MAX_NUMPROCS "")
set(KokkosCore_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(KokkosCore_ENABLE_Pthread OFF)
set(KokkosCore_ENABLE_CUDA OFF)
set(KokkosCore_ENABLE_HWLOC OFF)
set(KokkosCore_ENABLE_DLlib ON)

# Exported cache variables
set(Kokkos_CXX_STANDARD "")
set(Kokkos_ENABLE_THREADS "OFF")
set(Kokkos_ENABLE_OPENMP "OFF")
set(Kokkos_ENABLE_OPENACC "OFF")
set(Kokkos_ENABLE_OPENMPTARGET "OFF")
set(Kokkos_ENABLE_CUDA "OFF")
set(Kokkos_ENABLE_SERIAL "ON")
set(Kokkos_ENABLE_HPX "OFF")
set(Kokkos_ENABLE_HIP "OFF")
set(Kokkos_ENABLE_SYCL "OFF")
set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "OFF")
set(Kokkos_ENABLE_CUDA_UVM "OFF")
set(Kokkos_ENABLE_CUDA_LDG_INTRINSIC "OFF")
set(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "OFF")
set(Kokkos_ENABLE_DEPRECATED_CODE_3 "OFF")
set(Kokkos_ENABLE_DEPRECATED_CODE_4 "ON")
set(Kokkos_ENABLE_DEPRECATION_WARNINGS "ON")
set(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE "OFF")
set(Kokkos_ENABLE_HPX_ASYNC_DISPATCH "OFF")
set(Kokkos_ENABLE_BENCHMARKS "OFF")
set(Kokkos_ENABLE_DEBUG "OFF")
set(Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK "OFF")
set(Kokkos_ENABLE_LARGE_MEM_TESTS "OFF")
set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK "OFF")
set(Kokkos_ENABLE_COMPILER_WARNINGS "OFF")
set(Kokkos_ENABLE_PROFILING_LOAD_PRINT "OFF")
set(Kokkos_ENABLE_TUNING "OFF")
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION "OFF")
set(Kokkos_ENABLE_LAUNCH_COMPILER "ON")
set(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE "OFF")
set(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS "OFF")
set(Kokkos_ENABLE_IMPL_DESUL_ATOMICS "ON")
set(Kokkos_ENABLE_DESUL_ATOMICS_EXTERNAL "OFF")
set(Kokkos_ENABLE_IMPL_MDSPAN "OFF")
set(Kokkos_ENABLE_MDSPAN_EXTERNAL "OFF")
set(Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN "OFF")
set(Kokkos_ENABLE_CUDA_LAMBDA "OFF")
set(Kokkos_ENABLE_COMPLEX_ALIGN "OFF")
set(Kokkos_ENABLE_HEADER_SELF_CONTAINMENT_TESTS "OFF")
set(Kokkos_ENABLE_CUDA_CONSTEXPR "OFF")
set(Kokkos_ENABLE_UNSUPPORTED_ARCHS "OFF")
set(Kokkos_ENABLE_HWLOC "Off")
set(Kokkos_HWLOC_DIR "")
set(Kokkos_ENABLE_LIBNUMA "Off")
set(Kokkos_LIBNUMA_DIR "")
set(Kokkos_ENABLE_MEMKIND "Off")
set(Kokkos_MEMKIND_DIR "")
set(Kokkos_ENABLE_CUDA "OFF")
set(Kokkos_CUDA_DIR "")
set(Kokkos_ENABLE_LIBRT "Off")
set(Kokkos_LIBRT_DIR "")
set(Kokkos_ENABLE_ROCM "OFF")
set(Kokkos_ROCM_DIR "")
set(Kokkos_ENABLE_LIBDL "ON")
set(Kokkos_LIBDL_DIR "")
set(Kokkos_ENABLE_HPX "OFF")
set(Kokkos_HPX_DIR "")
set(Kokkos_ENABLE_THREADS "OFF")
set(Kokkos_THREADS_DIR "")
set(Kokkos_ENABLE_LIBQUADMATH "OFF")
set(Kokkos_LIBQUADMATH_DIR "")

# Include configuration of dependent packages
if (NOT TARGET DLlib::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/DLlib/DLlibConfig.cmake")
endif()

# Import KokkosCore targets
include("${CMAKE_CURRENT_LIST_DIR}/KokkosCoreTargets.cmake")

# Standard TriBITS-compliant external package variables
set(KokkosCore_IS_TRIBITS_COMPLIANT TRUE)
set(KokkosCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(KokkosCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(KokkosCore_EXPORTED_PACKAGE_LIBS_NAMES "kokkoscore")

foreach(libname IN LISTS KokkosCore_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE KokkosCore::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'KokkosCore::${libname}', or better yet,"
      " 'KokkosCore::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'KokkosCore'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'KokkosCore_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
