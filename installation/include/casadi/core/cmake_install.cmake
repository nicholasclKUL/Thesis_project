# Install script for directory: /home/nicholas/casadi/casadi/core

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6"
         RPATH "/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/nicholas/casadi/build/lib/libcasadi.so.3.6")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6"
         OLD_RPATH "::::::::::::::"
         NEW_RPATH "/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so.3.6")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so"
         RPATH "/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/nicholas/casadi/build/lib/libcasadi.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so"
         OLD_RPATH "::::::::::::::"
         NEW_RPATH "/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcasadi.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/casadi/core/casadi_limits.hpp;/usr/local/include/casadi/core/casadi_types.hpp;/usr/local/include/casadi/core/casadi_common.hpp;/usr/local/include/casadi/core/casadi_logger.hpp;/usr/local/include/casadi/core/casadi_interrupt.hpp;/usr/local/include/casadi/core/exception.hpp;/usr/local/include/casadi/core/calculus.hpp;/usr/local/include/casadi/core/global_options.hpp;/usr/local/include/casadi/core/casadi_meta.hpp;/usr/local/include/casadi/core/printable.hpp;/usr/local/include/casadi/core/shared_object.hpp;/usr/local/include/casadi/core/generic_type.hpp;/usr/local/include/casadi/core/options.hpp;/usr/local/include/casadi/core/casadi_misc.hpp;/usr/local/include/casadi/core/timing.hpp;/usr/local/include/casadi/core/polynomial.hpp;/usr/local/include/casadi/core/generic_expression.hpp;/usr/local/include/casadi/core/generic_matrix.hpp;/usr/local/include/casadi/core/matrix_fwd.hpp;/usr/local/include/casadi/core/matrix_decl.hpp;/usr/local/include/casadi/core/sx_fwd.hpp;/usr/local/include/casadi/core/sx.hpp;/usr/local/include/casadi/core/dm_fwd.hpp;/usr/local/include/casadi/core/dm.hpp;/usr/local/include/casadi/core/im_fwd.hpp;/usr/local/include/casadi/core/im.hpp;/usr/local/include/casadi/core/sparsity_interface.hpp;/usr/local/include/casadi/core/sparsity.hpp;/usr/local/include/casadi/core/slice.hpp;/usr/local/include/casadi/core/submatrix.hpp;/usr/local/include/casadi/core/nonzeros.hpp;/usr/local/include/casadi/core/sx_elem.hpp;/usr/local/include/casadi/core/sx.hpp;/usr/local/include/casadi/core/mx.hpp;/usr/local/include/casadi/core/function.hpp;/usr/local/include/casadi/core/callback.hpp;/usr/local/include/casadi/core/external.hpp;/usr/local/include/casadi/core/linsol.hpp;/usr/local/include/casadi/core/rootfinder.hpp;/usr/local/include/casadi/core/integrator.hpp;/usr/local/include/casadi/core/nlpsol.hpp;/usr/local/include/casadi/core/conic.hpp;/usr/local/include/casadi/core/dple.hpp;/usr/local/include/casadi/core/interpolant.hpp;/usr/local/include/casadi/core/expm.hpp;/usr/local/include/casadi/core/code_generator.hpp;/usr/local/include/casadi/core/importer.hpp;/usr/local/include/casadi/core/integration_tools.hpp;/usr/local/include/casadi/core/nlp_tools.hpp;/usr/local/include/casadi/core/nlp_builder.hpp;/usr/local/include/casadi/core/xml_node.hpp;/usr/local/include/casadi/core/xml_file.hpp;/usr/local/include/casadi/core/variable.hpp;/usr/local/include/casadi/core/dae_builder.hpp;/usr/local/include/casadi/core/optistack.hpp;/usr/local/include/casadi/core/serializer.hpp;/usr/local/include/casadi/core/serializing_stream.hpp;/usr/local/include/casadi/core/core.hpp;/usr/local/include/casadi/core/casadi_export.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/include/casadi/core" TYPE FILE FILES
    "/home/nicholas/casadi/casadi/core/casadi_limits.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_types.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_common.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_logger.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_interrupt.hpp"
    "/home/nicholas/casadi/casadi/core/exception.hpp"
    "/home/nicholas/casadi/casadi/core/calculus.hpp"
    "/home/nicholas/casadi/casadi/core/global_options.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_meta.hpp"
    "/home/nicholas/casadi/casadi/core/printable.hpp"
    "/home/nicholas/casadi/casadi/core/shared_object.hpp"
    "/home/nicholas/casadi/casadi/core/generic_type.hpp"
    "/home/nicholas/casadi/casadi/core/options.hpp"
    "/home/nicholas/casadi/casadi/core/casadi_misc.hpp"
    "/home/nicholas/casadi/casadi/core/timing.hpp"
    "/home/nicholas/casadi/casadi/core/polynomial.hpp"
    "/home/nicholas/casadi/casadi/core/generic_expression.hpp"
    "/home/nicholas/casadi/casadi/core/generic_matrix.hpp"
    "/home/nicholas/casadi/casadi/core/matrix_fwd.hpp"
    "/home/nicholas/casadi/casadi/core/matrix_decl.hpp"
    "/home/nicholas/casadi/casadi/core/sx_fwd.hpp"
    "/home/nicholas/casadi/casadi/core/sx.hpp"
    "/home/nicholas/casadi/casadi/core/dm_fwd.hpp"
    "/home/nicholas/casadi/casadi/core/dm.hpp"
    "/home/nicholas/casadi/casadi/core/im_fwd.hpp"
    "/home/nicholas/casadi/casadi/core/im.hpp"
    "/home/nicholas/casadi/casadi/core/sparsity_interface.hpp"
    "/home/nicholas/casadi/casadi/core/sparsity.hpp"
    "/home/nicholas/casadi/casadi/core/slice.hpp"
    "/home/nicholas/casadi/casadi/core/submatrix.hpp"
    "/home/nicholas/casadi/casadi/core/nonzeros.hpp"
    "/home/nicholas/casadi/casadi/core/sx_elem.hpp"
    "/home/nicholas/casadi/casadi/core/sx.hpp"
    "/home/nicholas/casadi/casadi/core/mx.hpp"
    "/home/nicholas/casadi/casadi/core/function.hpp"
    "/home/nicholas/casadi/casadi/core/callback.hpp"
    "/home/nicholas/casadi/casadi/core/external.hpp"
    "/home/nicholas/casadi/casadi/core/linsol.hpp"
    "/home/nicholas/casadi/casadi/core/rootfinder.hpp"
    "/home/nicholas/casadi/casadi/core/integrator.hpp"
    "/home/nicholas/casadi/casadi/core/nlpsol.hpp"
    "/home/nicholas/casadi/casadi/core/conic.hpp"
    "/home/nicholas/casadi/casadi/core/dple.hpp"
    "/home/nicholas/casadi/casadi/core/interpolant.hpp"
    "/home/nicholas/casadi/casadi/core/expm.hpp"
    "/home/nicholas/casadi/casadi/core/code_generator.hpp"
    "/home/nicholas/casadi/casadi/core/importer.hpp"
    "/home/nicholas/casadi/casadi/core/integration_tools.hpp"
    "/home/nicholas/casadi/casadi/core/nlp_tools.hpp"
    "/home/nicholas/casadi/casadi/core/nlp_builder.hpp"
    "/home/nicholas/casadi/casadi/core/xml_node.hpp"
    "/home/nicholas/casadi/casadi/core/xml_file.hpp"
    "/home/nicholas/casadi/casadi/core/variable.hpp"
    "/home/nicholas/casadi/casadi/core/dae_builder.hpp"
    "/home/nicholas/casadi/casadi/core/optistack.hpp"
    "/home/nicholas/casadi/casadi/core/serializer.hpp"
    "/home/nicholas/casadi/casadi/core/serializing_stream.hpp"
    "/home/nicholas/casadi/casadi/core/core.hpp"
    "/home/nicholas/casadi/build/casadi/core/casadi_export.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/nicholas/casadi/build/casadi/core/runtime/cmake_install.cmake")

endif()

