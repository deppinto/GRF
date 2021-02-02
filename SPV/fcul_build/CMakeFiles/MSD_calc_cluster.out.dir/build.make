# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/diogo/MEGA/cenas/CellGPU/SPV

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build

# Include any dependencies generated for this target.
include CMakeFiles/MSD_calc_cluster.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MSD_calc_cluster.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MSD_calc_cluster.out.dir/flags.make

CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o: CMakeFiles/MSD_calc_cluster.out.dir/flags.make
CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o: ../MSD_calc_cluster.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o -c /home/diogo/MEGA/cenas/CellGPU/SPV/MSD_calc_cluster.cpp

CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/diogo/MEGA/cenas/CellGPU/SPV/MSD_calc_cluster.cpp > CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.i

CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/diogo/MEGA/cenas/CellGPU/SPV/MSD_calc_cluster.cpp -o CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.s

# Object files for target MSD_calc_cluster.out
MSD_calc_cluster_out_OBJECTS = \
"CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o"

# External object files for target MSD_calc_cluster.out
MSD_calc_cluster_out_EXTERNAL_OBJECTS =

MSD_calc_cluster.out: CMakeFiles/MSD_calc_cluster.out.dir/MSD_calc_cluster.cpp.o
MSD_calc_cluster.out: CMakeFiles/MSD_calc_cluster.out.dir/build.make
MSD_calc_cluster.out: /usr/local/cuda-10.2/lib64/libcudart_static.a
MSD_calc_cluster.out: /usr/lib/x86_64-linux-gnu/librt.so
MSD_calc_cluster.out: src/models/libmodel.a
MSD_calc_cluster.out: src/models/libmodelGPU.a
MSD_calc_cluster.out: src/updaters/libupdaters.a
MSD_calc_cluster.out: src/updaters/libupdatersGPU.a
MSD_calc_cluster.out: src/analysis/libanalysis.a
MSD_calc_cluster.out: src/databases/libdatabase.a
MSD_calc_cluster.out: src/simulation/libsimulation.a
MSD_calc_cluster.out: src/utility/libutility.a
MSD_calc_cluster.out: src/utility/libutilityGPU.a
MSD_calc_cluster.out: /usr/local/lib/libCGAL.so.13.0.3
MSD_calc_cluster.out: /usr/local/cuda-10.2/lib64/libcudart_static.a
MSD_calc_cluster.out: /usr/lib/x86_64-linux-gnu/librt.so
MSD_calc_cluster.out: /usr/local/lib/libmpfr.so
MSD_calc_cluster.out: /usr/local/lib/libgmp.so
MSD_calc_cluster.out: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
MSD_calc_cluster.out: /usr/lib/x86_64-linux-gnu/libpthread.so
MSD_calc_cluster.out: CMakeFiles/MSD_calc_cluster.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MSD_calc_cluster.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MSD_calc_cluster.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MSD_calc_cluster.out.dir/build: MSD_calc_cluster.out

.PHONY : CMakeFiles/MSD_calc_cluster.out.dir/build

CMakeFiles/MSD_calc_cluster.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MSD_calc_cluster.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MSD_calc_cluster.out.dir/clean

CMakeFiles/MSD_calc_cluster.out.dir/depend:
	cd /home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/diogo/MEGA/cenas/CellGPU/SPV /home/diogo/MEGA/cenas/CellGPU/SPV /home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build /home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build /home/diogo/MEGA/cenas/CellGPU/SPV/fcul_build/CMakeFiles/MSD_calc_cluster.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MSD_calc_cluster.out.dir/depend
