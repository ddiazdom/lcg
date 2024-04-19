# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = /Users/ddiaz/Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Users/ddiaz/Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ddiaz/CLionProjects/lcg-2023

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/reverse_bwt.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/reverse_bwt.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/reverse_bwt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/reverse_bwt.dir/flags.make

CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o: CMakeFiles/reverse_bwt.dir/flags.make
CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o: ../scripts/reverse_bwt.cpp
CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o: CMakeFiles/reverse_bwt.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o -MF CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o.d -o CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o -c /Users/ddiaz/CLionProjects/lcg-2023/scripts/reverse_bwt.cpp

CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ddiaz/CLionProjects/lcg-2023/scripts/reverse_bwt.cpp > CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.i

CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ddiaz/CLionProjects/lcg-2023/scripts/reverse_bwt.cpp -o CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.s

# Object files for target reverse_bwt
reverse_bwt_OBJECTS = \
"CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o"

# External object files for target reverse_bwt
reverse_bwt_EXTERNAL_OBJECTS =

reverse_bwt: CMakeFiles/reverse_bwt.dir/scripts/reverse_bwt.cpp.o
reverse_bwt: CMakeFiles/reverse_bwt.dir/build.make
reverse_bwt: libgrlbwt.a
reverse_bwt: external/bioparsers/libbiopar.a
reverse_bwt: external/cdt/libcdt.a
reverse_bwt: /Users/ddiaz/lib/libsdsl.a
reverse_bwt: /Users/ddiaz/lib/libdivsufsort.a
reverse_bwt: /Users/ddiaz/lib/libdivsufsort64.a
reverse_bwt: CMakeFiles/reverse_bwt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable reverse_bwt"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/reverse_bwt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/reverse_bwt.dir/build: reverse_bwt
.PHONY : CMakeFiles/reverse_bwt.dir/build

CMakeFiles/reverse_bwt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/reverse_bwt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/reverse_bwt.dir/clean

CMakeFiles/reverse_bwt.dir/depend:
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ddiaz/CLionProjects/lcg-2023 /Users/ddiaz/CLionProjects/lcg-2023 /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles/reverse_bwt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/reverse_bwt.dir/depend
