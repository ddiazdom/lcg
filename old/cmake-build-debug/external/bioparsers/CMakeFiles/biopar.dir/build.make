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
include external/bioparsers/CMakeFiles/biopar.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/bioparsers/CMakeFiles/biopar.dir/compiler_depend.make

# Include the progress variables for this target.
include external/bioparsers/CMakeFiles/biopar.dir/progress.make

# Include the compile flags for this target's objects.
include external/bioparsers/CMakeFiles/biopar.dir/flags.make

external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o: external/bioparsers/CMakeFiles/biopar.dir/flags.make
external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o: ../external/bioparsers/lib/dna_string.cpp
external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o: external/bioparsers/CMakeFiles/biopar.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o -MF CMakeFiles/biopar.dir/lib/dna_string.cpp.o.d -o CMakeFiles/biopar.dir/lib/dna_string.cpp.o -c /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/dna_string.cpp

external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/biopar.dir/lib/dna_string.cpp.i"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/dna_string.cpp > CMakeFiles/biopar.dir/lib/dna_string.cpp.i

external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/biopar.dir/lib/dna_string.cpp.s"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/dna_string.cpp -o CMakeFiles/biopar.dir/lib/dna_string.cpp.s

external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o: external/bioparsers/CMakeFiles/biopar.dir/flags.make
external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o: ../external/bioparsers/lib/fastx_handler.cpp
external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o: external/bioparsers/CMakeFiles/biopar.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o -MF CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o.d -o CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o -c /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/fastx_handler.cpp

external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/biopar.dir/lib/fastx_handler.cpp.i"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/fastx_handler.cpp > CMakeFiles/biopar.dir/lib/fastx_handler.cpp.i

external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/biopar.dir/lib/fastx_handler.cpp.s"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers/lib/fastx_handler.cpp -o CMakeFiles/biopar.dir/lib/fastx_handler.cpp.s

# Object files for target biopar
biopar_OBJECTS = \
"CMakeFiles/biopar.dir/lib/dna_string.cpp.o" \
"CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o"

# External object files for target biopar
biopar_EXTERNAL_OBJECTS =

external/bioparsers/libbiopar.a: external/bioparsers/CMakeFiles/biopar.dir/lib/dna_string.cpp.o
external/bioparsers/libbiopar.a: external/bioparsers/CMakeFiles/biopar.dir/lib/fastx_handler.cpp.o
external/bioparsers/libbiopar.a: external/bioparsers/CMakeFiles/biopar.dir/build.make
external/bioparsers/libbiopar.a: external/bioparsers/CMakeFiles/biopar.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libbiopar.a"
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && $(CMAKE_COMMAND) -P CMakeFiles/biopar.dir/cmake_clean_target.cmake
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/biopar.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/bioparsers/CMakeFiles/biopar.dir/build: external/bioparsers/libbiopar.a
.PHONY : external/bioparsers/CMakeFiles/biopar.dir/build

external/bioparsers/CMakeFiles/biopar.dir/clean:
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers && $(CMAKE_COMMAND) -P CMakeFiles/biopar.dir/cmake_clean.cmake
.PHONY : external/bioparsers/CMakeFiles/biopar.dir/clean

external/bioparsers/CMakeFiles/biopar.dir/depend:
	cd /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ddiaz/CLionProjects/lcg-2023 /Users/ddiaz/CLionProjects/lcg-2023/external/bioparsers /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers /Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/external/bioparsers/CMakeFiles/biopar.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/bioparsers/CMakeFiles/biopar.dir/depend

