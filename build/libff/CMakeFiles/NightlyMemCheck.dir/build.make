# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /opt/cmake-3.14.1-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.14.1-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/newton/github/test/FFT-Multiexp-AMD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/newton/github/test/FFT-Multiexp-AMD/build

# Utility rule file for NightlyMemCheck.

# Include the progress variables for this target.
include libff/CMakeFiles/NightlyMemCheck.dir/progress.make

libff/CMakeFiles/NightlyMemCheck:
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && /opt/cmake-3.14.1-Linux-x86_64/bin/ctest -D NightlyMemCheck

NightlyMemCheck: libff/CMakeFiles/NightlyMemCheck
NightlyMemCheck: libff/CMakeFiles/NightlyMemCheck.dir/build.make

.PHONY : NightlyMemCheck

# Rule to build all files generated by this target.
libff/CMakeFiles/NightlyMemCheck.dir/build: NightlyMemCheck

.PHONY : libff/CMakeFiles/NightlyMemCheck.dir/build

libff/CMakeFiles/NightlyMemCheck.dir/clean:
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && $(CMAKE_COMMAND) -P CMakeFiles/NightlyMemCheck.dir/cmake_clean.cmake
.PHONY : libff/CMakeFiles/NightlyMemCheck.dir/clean

libff/CMakeFiles/NightlyMemCheck.dir/depend:
	cd /home/newton/github/test/FFT-Multiexp-AMD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/newton/github/test/FFT-Multiexp-AMD /home/newton/github/test/FFT-Multiexp-AMD/libff /home/newton/github/test/FFT-Multiexp-AMD/build /home/newton/github/test/FFT-Multiexp-AMD/build/libff /home/newton/github/test/FFT-Multiexp-AMD/build/libff/CMakeFiles/NightlyMemCheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libff/CMakeFiles/NightlyMemCheck.dir/depend

