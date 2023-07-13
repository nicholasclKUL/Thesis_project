# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nicholas/casadi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nicholas/casadi/build

# Utility rule file for libs.

# Include any custom commands dependencies for this target.
include casadi/CMakeFiles/libs.dir/compiler_depend.make

# Include the progress variables for this target.
include casadi/CMakeFiles/libs.dir/progress.make

casadi/CMakeFiles/libs: lib/libcasadi.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_sundials_common.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_integrator_cvodes.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_integrator_idas.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_rootfinder_kinsol.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_csparse.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_csparsecholesky.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_xmlfile_tinyxml.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_conic_nlpsol.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_conic_qrqp.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_nlpsol_qrsqp.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_importer_shell.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_integrator_rk.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_integrator_collocation.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_interpolant_linear.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_interpolant_bspline.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_symbolicqr.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_qr.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_ldl.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_tridiag.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_linsol_lsqr.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_nlpsol_sqpmethod.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_nlpsol_scpgen.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_rootfinder_newton.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_rootfinder_fast_newton.so.3.6
casadi/CMakeFiles/libs: lib/libcasadi_rootfinder_nlpsol.so.3.6

libs: casadi/CMakeFiles/libs
libs: casadi/CMakeFiles/libs.dir/build.make
.PHONY : libs

# Rule to build all files generated by this target.
casadi/CMakeFiles/libs.dir/build: libs
.PHONY : casadi/CMakeFiles/libs.dir/build

casadi/CMakeFiles/libs.dir/clean:
	cd /home/nicholas/casadi/build/casadi && $(CMAKE_COMMAND) -P CMakeFiles/libs.dir/cmake_clean.cmake
.PHONY : casadi/CMakeFiles/libs.dir/clean

casadi/CMakeFiles/libs.dir/depend:
	cd /home/nicholas/casadi/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicholas/casadi /home/nicholas/casadi/casadi /home/nicholas/casadi/build /home/nicholas/casadi/build/casadi /home/nicholas/casadi/build/casadi/CMakeFiles/libs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : casadi/CMakeFiles/libs.dir/depend
