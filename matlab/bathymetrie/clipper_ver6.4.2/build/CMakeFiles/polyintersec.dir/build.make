# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build

# Include any dependencies generated for this target.
include CMakeFiles/polyintersec.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/polyintersec.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/polyintersec.dir/flags.make

CMakeFiles/polyintersec.dir/src/clipper.cpp.o: CMakeFiles/polyintersec.dir/flags.make
CMakeFiles/polyintersec.dir/src/clipper.cpp.o: ../src/clipper.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/polyintersec.dir/src/clipper.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/polyintersec.dir/src/clipper.cpp.o -c /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/clipper.cpp

CMakeFiles/polyintersec.dir/src/clipper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polyintersec.dir/src/clipper.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/clipper.cpp > CMakeFiles/polyintersec.dir/src/clipper.cpp.i

CMakeFiles/polyintersec.dir/src/clipper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polyintersec.dir/src/clipper.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/clipper.cpp -o CMakeFiles/polyintersec.dir/src/clipper.cpp.s

CMakeFiles/polyintersec.dir/src/clipper.cpp.o.requires:
.PHONY : CMakeFiles/polyintersec.dir/src/clipper.cpp.o.requires

CMakeFiles/polyintersec.dir/src/clipper.cpp.o.provides: CMakeFiles/polyintersec.dir/src/clipper.cpp.o.requires
	$(MAKE) -f CMakeFiles/polyintersec.dir/build.make CMakeFiles/polyintersec.dir/src/clipper.cpp.o.provides.build
.PHONY : CMakeFiles/polyintersec.dir/src/clipper.cpp.o.provides

CMakeFiles/polyintersec.dir/src/clipper.cpp.o.provides.build: CMakeFiles/polyintersec.dir/src/clipper.cpp.o

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o: CMakeFiles/polyintersec.dir/flags.make
CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o: ../src/poly_intersection.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o -c /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/poly_intersection.cpp

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/poly_intersection.cpp > CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.i

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/src/poly_intersection.cpp -o CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.s

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.requires:
.PHONY : CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.requires

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.provides: CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.requires
	$(MAKE) -f CMakeFiles/polyintersec.dir/build.make CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.provides.build
.PHONY : CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.provides

CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.provides.build: CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o

# Object files for target polyintersec
polyintersec_OBJECTS = \
"CMakeFiles/polyintersec.dir/src/clipper.cpp.o" \
"CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o"

# External object files for target polyintersec
polyintersec_EXTERNAL_OBJECTS =

/localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec: CMakeFiles/polyintersec.dir/src/clipper.cpp.o
/localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec: CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o
/localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec: CMakeFiles/polyintersec.dir/build.make
/localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec: CMakeFiles/polyintersec.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable /localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polyintersec.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/polyintersec.dir/build: /localdata/gnautic/git/galileonautic/matlab/bathymetrie/polyintersec
.PHONY : CMakeFiles/polyintersec.dir/build

CMakeFiles/polyintersec.dir/requires: CMakeFiles/polyintersec.dir/src/clipper.cpp.o.requires
CMakeFiles/polyintersec.dir/requires: CMakeFiles/polyintersec.dir/src/poly_intersection.cpp.o.requires
.PHONY : CMakeFiles/polyintersec.dir/requires

CMakeFiles/polyintersec.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/polyintersec.dir/cmake_clean.cmake
.PHONY : CMakeFiles/polyintersec.dir/clean

CMakeFiles/polyintersec.dir/depend:
	cd /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2 /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2 /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build /localdata/gnautic/git/galileonautic/matlab/bathymetrie/clipper_ver6.4.2/build/CMakeFiles/polyintersec.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/polyintersec.dir/depend

