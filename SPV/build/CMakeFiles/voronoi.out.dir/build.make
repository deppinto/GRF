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
CMAKE_SOURCE_DIR = /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build

# Include any dependencies generated for this target.
include CMakeFiles/voronoi.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/voronoi.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/voronoi.out.dir/flags.make

CMakeFiles/voronoi.out.dir/voronoi.cpp.o: CMakeFiles/voronoi.out.dir/flags.make
CMakeFiles/voronoi.out.dir/voronoi.cpp.o: ../voronoi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/voronoi.out.dir/voronoi.cpp.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/voronoi.out.dir/voronoi.cpp.o -c /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/voronoi.cpp

CMakeFiles/voronoi.out.dir/voronoi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/voronoi.out.dir/voronoi.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/voronoi.cpp > CMakeFiles/voronoi.out.dir/voronoi.cpp.i

CMakeFiles/voronoi.out.dir/voronoi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/voronoi.out.dir/voronoi.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/voronoi.cpp -o CMakeFiles/voronoi.out.dir/voronoi.cpp.s

CMakeFiles/voronoi.out.dir/voronoi.cpp.o.requires:

.PHONY : CMakeFiles/voronoi.out.dir/voronoi.cpp.o.requires

CMakeFiles/voronoi.out.dir/voronoi.cpp.o.provides: CMakeFiles/voronoi.out.dir/voronoi.cpp.o.requires
	$(MAKE) -f CMakeFiles/voronoi.out.dir/build.make CMakeFiles/voronoi.out.dir/voronoi.cpp.o.provides.build
.PHONY : CMakeFiles/voronoi.out.dir/voronoi.cpp.o.provides

CMakeFiles/voronoi.out.dir/voronoi.cpp.o.provides.build: CMakeFiles/voronoi.out.dir/voronoi.cpp.o


# Object files for target voronoi.out
voronoi_out_OBJECTS = \
"CMakeFiles/voronoi.out.dir/voronoi.cpp.o"

# External object files for target voronoi.out
voronoi_out_EXTERNAL_OBJECTS =

voronoi.out: CMakeFiles/voronoi.out.dir/voronoi.cpp.o
voronoi.out: CMakeFiles/voronoi.out.dir/build.make
voronoi.out: /usr/local/cuda/lib64/libcudart_static.a
voronoi.out: /usr/lib/x86_64-linux-gnu/librt.so
voronoi.out: src/models/libmodel.a
voronoi.out: src/models/libmodelGPU.a
voronoi.out: src/updaters/libupdaters.a
voronoi.out: src/updaters/libupdatersGPU.a
voronoi.out: src/analysis/libanalysis.a
voronoi.out: src/databases/libdatabase.a
voronoi.out: src/simulation/libsimulation.a
voronoi.out: src/utility/libutility.a
voronoi.out: src/utility/libutilityGPU.a
voronoi.out: /usr/lib/x86_64-linux-gnu/libCGAL.so.13.0.1
voronoi.out: /usr/local/cuda/lib64/libcudart_static.a
voronoi.out: /usr/lib/x86_64-linux-gnu/librt.so
voronoi.out: /usr/lib/x86_64-linux-gnu/libgmp.so
voronoi.out: /usr/lib/x86_64-linux-gnu/libgmpxx.so
voronoi.out: /usr/lib/x86_64-linux-gnu/libmpfr.so
voronoi.out: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
voronoi.out: /usr/lib/x86_64-linux-gnu/libpthread.so
voronoi.out: CMakeFiles/voronoi.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable voronoi.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/voronoi.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/voronoi.out.dir/build: voronoi.out

.PHONY : CMakeFiles/voronoi.out.dir/build

CMakeFiles/voronoi.out.dir/requires: CMakeFiles/voronoi.out.dir/voronoi.cpp.o.requires

.PHONY : CMakeFiles/voronoi.out.dir/requires

CMakeFiles/voronoi.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/voronoi.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/voronoi.out.dir/clean

CMakeFiles/voronoi.out.dir/depend:
	cd /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles/voronoi.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/voronoi.out.dir/depend
