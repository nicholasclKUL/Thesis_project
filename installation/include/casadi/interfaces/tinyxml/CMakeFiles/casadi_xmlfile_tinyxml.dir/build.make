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

# Include any dependencies generated for this target.
include casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/compiler_depend.make

# Include the progress variables for this target.
include casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/progress.make

# Include the compile flags for this target's objects.
include casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/flags.make

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/flags.make
casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o: ../casadi/interfaces/tinyxml/tinyxml_interface.cpp
casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicholas/casadi/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o -MF CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o.d -o CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o -c /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface.cpp

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.i"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface.cpp > CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.i

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.s"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface.cpp -o CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.s

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/flags.make
casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o: ../casadi/interfaces/tinyxml/tinyxml_interface_meta.cpp
casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicholas/casadi/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o -MF CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o.d -o CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o -c /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface_meta.cpp

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.i"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface_meta.cpp > CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.i

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.s"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicholas/casadi/casadi/interfaces/tinyxml/tinyxml_interface_meta.cpp -o CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.s

# Object files for target casadi_xmlfile_tinyxml
casadi_xmlfile_tinyxml_OBJECTS = \
"CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o" \
"CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o"

# External object files for target casadi_xmlfile_tinyxml
casadi_xmlfile_tinyxml_EXTERNAL_OBJECTS =

lib/libcasadi_xmlfile_tinyxml.so.3.6: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface.cpp.o
lib/libcasadi_xmlfile_tinyxml.so.3.6: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/tinyxml_interface_meta.cpp.o
lib/libcasadi_xmlfile_tinyxml.so.3.6: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/build.make
lib/libcasadi_xmlfile_tinyxml.so.3.6: lib/libcasadi.so.3.6
lib/libcasadi_xmlfile_tinyxml.so.3.6: lib/libcasadi_tinyxml.a
lib/libcasadi_xmlfile_tinyxml.so.3.6: casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nicholas/casadi/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library ../../../lib/libcasadi_xmlfile_tinyxml.so"
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/casadi_xmlfile_tinyxml.dir/link.txt --verbose=$(VERBOSE)
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && $(CMAKE_COMMAND) -E cmake_symlink_library ../../../lib/libcasadi_xmlfile_tinyxml.so.3.6 ../../../lib/libcasadi_xmlfile_tinyxml.so.3.6 ../../../lib/libcasadi_xmlfile_tinyxml.so

lib/libcasadi_xmlfile_tinyxml.so: lib/libcasadi_xmlfile_tinyxml.so.3.6
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libcasadi_xmlfile_tinyxml.so

# Rule to build all files generated by this target.
casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/build: lib/libcasadi_xmlfile_tinyxml.so
.PHONY : casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/build

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/clean:
	cd /home/nicholas/casadi/build/casadi/interfaces/tinyxml && $(CMAKE_COMMAND) -P CMakeFiles/casadi_xmlfile_tinyxml.dir/cmake_clean.cmake
.PHONY : casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/clean

casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/depend:
	cd /home/nicholas/casadi/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicholas/casadi /home/nicholas/casadi/casadi/interfaces/tinyxml /home/nicholas/casadi/build /home/nicholas/casadi/build/casadi/interfaces/tinyxml /home/nicholas/casadi/build/casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : casadi/interfaces/tinyxml/CMakeFiles/casadi_xmlfile_tinyxml.dir/depend

