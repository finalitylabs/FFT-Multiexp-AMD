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

# Include any dependencies generated for this target.
include libff/CMakeFiles/multiexp_profile.dir/depend.make

# Include the progress variables for this target.
include libff/CMakeFiles/multiexp_profile.dir/progress.make

# Include the compile flags for this target's objects.
include libff/CMakeFiles/multiexp_profile.dir/flags.make

libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o: libff/CMakeFiles/multiexp_profile.dir/flags.make
libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o: ../libff/algebra/scalar_multiplication/multiexp_profile.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/newton/github/test/FFT-Multiexp-AMD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o"
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o -c /home/newton/github/test/FFT-Multiexp-AMD/libff/algebra/scalar_multiplication/multiexp_profile.cpp

libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.i"
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/newton/github/test/FFT-Multiexp-AMD/libff/algebra/scalar_multiplication/multiexp_profile.cpp > CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.i

libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.s"
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/newton/github/test/FFT-Multiexp-AMD/libff/algebra/scalar_multiplication/multiexp_profile.cpp -o CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.s

# Object files for target multiexp_profile
multiexp_profile_OBJECTS = \
"CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o"

# External object files for target multiexp_profile
multiexp_profile_EXTERNAL_OBJECTS =

libff/multiexp_profile: libff/CMakeFiles/multiexp_profile.dir/algebra/scalar_multiplication/multiexp_profile.cpp.o
libff/multiexp_profile: libff/CMakeFiles/multiexp_profile.dir/build.make
libff/multiexp_profile: libff/libff.a
libff/multiexp_profile: /usr/lib/x86_64-linux-gnu/libgmp.so
libff/multiexp_profile: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libff/multiexp_profile: depends/libzm.a
libff/multiexp_profile: libff/CMakeFiles/multiexp_profile.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/newton/github/test/FFT-Multiexp-AMD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable multiexp_profile"
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/multiexp_profile.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libff/CMakeFiles/multiexp_profile.dir/build: libff/multiexp_profile

.PHONY : libff/CMakeFiles/multiexp_profile.dir/build

libff/CMakeFiles/multiexp_profile.dir/clean:
	cd /home/newton/github/test/FFT-Multiexp-AMD/build/libff && $(CMAKE_COMMAND) -P CMakeFiles/multiexp_profile.dir/cmake_clean.cmake
.PHONY : libff/CMakeFiles/multiexp_profile.dir/clean

libff/CMakeFiles/multiexp_profile.dir/depend:
	cd /home/newton/github/test/FFT-Multiexp-AMD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/newton/github/test/FFT-Multiexp-AMD /home/newton/github/test/FFT-Multiexp-AMD/libff /home/newton/github/test/FFT-Multiexp-AMD/build /home/newton/github/test/FFT-Multiexp-AMD/build/libff /home/newton/github/test/FFT-Multiexp-AMD/build/libff/CMakeFiles/multiexp_profile.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libff/CMakeFiles/multiexp_profile.dir/depend

