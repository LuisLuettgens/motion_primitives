# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/luis/Backup2018Seafile/galileonautic

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/luis/Backup2018Seafile/galileonautic/build

# Include any dependencies generated for this target.
include external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/depend.make

# Include the progress variables for this target.
include external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/progress.make

# Include the compile flags for this target's objects.
include external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o: ../external/transworhp/src/toolbase/tool.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/tool.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/tool.cpp > CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/tool.cpp -o CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o: ../external/transworhp/src/toolbase/toolstatus.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolstatus.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolstatus.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolstatus.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o: ../external/transworhp/src/toolbase/toolmenu.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenu.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenu.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenu.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o: ../external/transworhp/src/toolbase/toolrect.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolrect.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolrect.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolrect.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o: ../external/transworhp/src/toolbase/toolmenuentry.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenuentry.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenuentry.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenuentry.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o: ../external/transworhp/src/toolbase/toolmenukey.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenukey.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenukey.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenukey.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o


external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/flags.make
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o: ../external/transworhp/src/toolbase/toolwindow.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/luis/Backup2018Seafile/galileonautic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o -c /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolwindow.cpp

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.i"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolwindow.cpp > CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.i

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.s"
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolwindow.cpp -o CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.s

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.requires:

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.provides: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.requires
	$(MAKE) -f external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.provides.build
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.provides

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.provides.build: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o


TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o
TransWORHP_toolbase: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build.make

.PHONY : TransWORHP_toolbase

# Rule to build all files generated by this target.
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build: TransWORHP_toolbase

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/build

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/tool.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolstatus.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenu.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolrect.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenuentry.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolmenukey.cpp.o.requires
external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires: external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/toolwindow.cpp.o.requires

.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/requires

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/clean:
	cd /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase && $(CMAKE_COMMAND) -P CMakeFiles/TransWORHP_toolbase.dir/cmake_clean.cmake
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/clean

external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/depend:
	cd /home/luis/Backup2018Seafile/galileonautic/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/luis/Backup2018Seafile/galileonautic /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase /home/luis/Backup2018Seafile/galileonautic/build /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase /home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/transworhp/src/toolbase/CMakeFiles/TransWORHP_toolbase.dir/depend

