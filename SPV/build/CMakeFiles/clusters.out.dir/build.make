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
include CMakeFiles/clusters.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/clusters.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/clusters.out.dir/flags.make

CMakeFiles/clusters.out.dir/clusters.cpp.o: CMakeFiles/clusters.out.dir/flags.make
CMakeFiles/clusters.out.dir/clusters.cpp.o: ../clusters.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/clusters.out.dir/clusters.cpp.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/clusters.out.dir/clusters.cpp.o -c /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/clusters.cpp

CMakeFiles/clusters.out.dir/clusters.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/clusters.out.dir/clusters.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/clusters.cpp > CMakeFiles/clusters.out.dir/clusters.cpp.i

CMakeFiles/clusters.out.dir/clusters.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/clusters.out.dir/clusters.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/clusters.cpp -o CMakeFiles/clusters.out.dir/clusters.cpp.s

CMakeFiles/clusters.out.dir/clusters.cpp.o.requires:

.PHONY : CMakeFiles/clusters.out.dir/clusters.cpp.o.requires

CMakeFiles/clusters.out.dir/clusters.cpp.o.provides: CMakeFiles/clusters.out.dir/clusters.cpp.o.requires
	$(MAKE) -f CMakeFiles/clusters.out.dir/build.make CMakeFiles/clusters.out.dir/clusters.cpp.o.provides.build
.PHONY : CMakeFiles/clusters.out.dir/clusters.cpp.o.provides

CMakeFiles/clusters.out.dir/clusters.cpp.o.provides.build: CMakeFiles/clusters.out.dir/clusters.cpp.o


# Object files for target clusters.out
clusters_out_OBJECTS = \
"CMakeFiles/clusters.out.dir/clusters.cpp.o"

# External object files for target clusters.out
clusters_out_EXTERNAL_OBJECTS =

clusters.out: CMakeFiles/clusters.out.dir/clusters.cpp.o
clusters.out: CMakeFiles/clusters.out.dir/build.make
clusters.out: /usr/local/cuda/lib64/libcudart_static.a
clusters.out: /usr/lib/x86_64-linux-gnu/librt.so
clusters.out: src/models/libmodel.a
clusters.out: src/models/libmodelGPU.a
clusters.out: src/updaters/libupdaters.a
clusters.out: src/updaters/libupdatersGPU.a
clusters.out: src/analysis/libanalysis.a
clusters.out: src/databases/libdatabase.a
clusters.out: src/simulation/libsimulation.a
clusters.out: src/utility/libutility.a
clusters.out: src/utility/libutilityGPU.a
clusters.out: /usr/lib/x86_64-linux-gnu/libCGAL.so.13.0.1
clusters.out: /usr/local/cuda/lib64/libcudart_static.a
clusters.out: /usr/lib/x86_64-linux-gnu/librt.so
clusters.out: /usr/lib/x86_64-linux-gnu/libgmp.so
clusters.out: /usr/lib/x86_64-linux-gnu/libgmpxx.so
clusters.out: /usr/lib/x86_64-linux-gnu/libmpfr.so
clusters.out: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
clusters.out: /usr/lib/x86_64-linux-gnu/libpthread.so
clusters.out: CMakeFiles/clusters.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable clusters.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/clusters.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/clusters.out.dir/build: clusters.out

.PHONY : CMakeFiles/clusters.out.dir/build

CMakeFiles/clusters.out.dir/requires: CMakeFiles/clusters.out.dir/clusters.cpp.o.requires

.PHONY : CMakeFiles/clusters.out.dir/requires

CMakeFiles/clusters.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/clusters.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/clusters.out.dir/clean

CMakeFiles/clusters.out.dir/depend:
	cd /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles/clusters.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/clusters.out.dir/depend

