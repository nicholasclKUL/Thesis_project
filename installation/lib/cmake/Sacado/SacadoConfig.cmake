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
# CMake variable for use by Trilinos/Sacado clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Sacado requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Sacado build
## ---------------------------------------------------------------------------

set(Sacado_CXX_COMPILER "/usr/bin/g++")

set(Sacado_C_COMPILER "/usr/bin/gcc")

set(Sacado_Fortran_COMPILER "/usr/bin/gfortran")
# Deprecated!
set(Sacado_FORTRAN_COMPILER "/usr/bin/gfortran") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Sacado build
## ---------------------------------------------------------------------------

## Give the build type
set(Sacado_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Sacado_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Sacado_C_FLAGS [[ -pedantic -Wall -Wno-long-long -std=c99 -O3 -DNDEBUG]])

set(Sacado_Fortran_FLAGS [[ -O3]])
# Deprecated
set(Sacado_FORTRAN_FLAGS [[ -O3]])

## Extra link flags (e.g., specification of fortran libraries)
set(Sacado_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Sacado_SHARED_LIB_RPATH_COMMAND "")
set(Sacado_BUILD_SHARED_LIBS "FALSE")

set(Sacado_LINKER /usr/bin/ld)
set(Sacado_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Sacado_INSTALL_DIR "/home/nicholas/trilinos")

## List of package libraries
set(Sacado_LIBRARIES Sacado::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Sacado_MPI_LIBRARIES "")
set(Sacado_MPI_LIBRARY_DIRS "")
set(Sacado_MPI_INCLUDE_DIRS "")
set(Sacado_MPI_EXEC "")
set(Sacado_MPI_EXEC_MAX_NUMPROCS "")
set(Sacado_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Sacado_ENABLE_KokkosCore ON)
set(Sacado_ENABLE_TeuchosCore ON)
set(Sacado_ENABLE_TeuchosNumerics OFF)
set(Sacado_ENABLE_TeuchosComm ON)
set(Sacado_ENABLE_TeuchosKokkosComm ON)

# Exported cache variables
set(Sacado_ENABLE_COMPLEX "FALSE")
set(HAVE_SACADO_COMPLEX "OFF")
set(Sacado_NEW_FAD_DESIGN_IS_DEFAULT "ON")
set(SACADO_NEW_FAD_DESIGN_IS_DEFAULT "ON")
set(Sacado_SFAD_INIT_DEFAULT_CONSTRUCTOR "OFF")
set(SACADO_SFAD_INIT_DEFAULT_CONSTRUCTOR "OFF")
set(Sacado_ENABLE_HIERARCHICAL "FALSE")
set(SACADO_VIEW_CUDA_HIERARCHICAL "OFF")
set(Sacado_ENABLE_HIERARCHICAL_DFAD "FALSE")
set(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD "OFF")
set(Sacado_ENABLE_MEMORY_POOL "FALSE")
set(SACADO_KOKKOS_USE_MEMORY_POOL "OFF")
set(Sacado_ENABLE_DEBUG "OFF")
set(SACADO_DEBUG "OFF")
set(Sacado_ENABLE_UNINIT "OFF")
set(HAVE_SACADO_UNINIT "OFF")
set(Sacado_ENABLE_VIEW_SPEC "ON")
set(HAVE_SACADO_VIEW_SPEC "ON")
set(Sacado_ENABLE_GTest "ON")
set(HAVE_SACADO_GTEST "ON")
set(Sacado_ENABLE_ADOLC "")
set(HAVE_ADOLC "OFF")
set(Sacado_ENABLE_ADIC "")
set(HAVE_ADIC "OFF")
set(Sacado_ENABLE_RAD_NO_USING_STDCC "OFF")
set(RAD_NO_USING_STDCC "OFF")

# Include configuration of dependent packages
if (NOT TARGET KokkosCore::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../KokkosCore/KokkosCoreConfig.cmake")
endif()
if (NOT TARGET TeuchosCore::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosCore/TeuchosCoreConfig.cmake")
endif()
if (NOT TARGET TeuchosComm::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosComm/TeuchosCommConfig.cmake")
endif()
if (NOT TARGET TeuchosKokkosComm::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosKokkosComm/TeuchosKokkosCommConfig.cmake")
endif()

# Import Sacado targets
include("${CMAKE_CURRENT_LIST_DIR}/SacadoTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Sacado_IS_TRIBITS_COMPLIANT TRUE)
set(Sacado_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Sacado_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Sacado_EXPORTED_PACKAGE_LIBS_NAMES "sacado")

foreach(libname IN LISTS Sacado_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Sacado::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Sacado::${libname}', or better yet,"
      " 'Sacado::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Sacado'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Sacado_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
