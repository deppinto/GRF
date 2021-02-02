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
include CMakeFiles/percT.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/percT.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/percT.out.dir/flags.make

CMakeFiles/percT.out.dir/percT.cpp.o: CMakeFiles/percT.out.dir/flags.make
CMakeFiles/percT.out.dir/percT.cpp.o: ../percT.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/percT.out.dir/percT.cpp.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/percT.out.dir/percT.cpp.o -c /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/percT.cpp

CMakeFiles/percT.out.dir/percT.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/percT.out.dir/percT.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/percT.cpp > CMakeFiles/percT.out.dir/percT.cpp.i

CMakeFiles/percT.out.dir/percT.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/percT.out.dir/percT.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/percT.cpp -o CMakeFiles/percT.out.dir/percT.cpp.s

CMakeFiles/percT.out.dir/percT.cpp.o.requires:

.PHONY : CMakeFiles/percT.out.dir/percT.cpp.o.requires

CMakeFiles/percT.out.dir/percT.cpp.o.provides: CMakeFiles/percT.out.dir/percT.cpp.o.requires
	$(MAKE) -f CMakeFiles/percT.out.dir/build.make CMakeFiles/percT.out.dir/percT.cpp.o.provides.build
.PHONY : CMakeFiles/percT.out.dir/percT.cpp.o.provides

CMakeFiles/percT.out.dir/percT.cpp.o.provides.build: CMakeFiles/percT.out.dir/percT.cpp.o


# Object files for target percT.out
percT_out_OBJECTS = \
"CMakeFiles/percT.out.dir/percT.cpp.o"

# External object files for target percT.out
percT_out_EXTERNAL_OBJECTS =

percT.out: CMakeFiles/percT.out.dir/percT.cpp.o
percT.out: CMakeFiles/percT.out.dir/build.make
percT.out: /usr/local/cuda/lib64/libcudart_static.a
percT.out: /usr/lib/x86_64-linux-gnu/librt.so
percT.out: src/models/libmodel.a
percT.out: src/models/libmodelGPU.a
percT.out: src/updaters/libupdaters.a
percT.out: src/updaters/libupdatersGPU.a
percT.out: src/analysis/libanalysis.a
percT.out: src/databases/libdatabase.a
percT.out: src/simulation/libsimulation.a
percT.out: src/utility/libutility.a
percT.out: src/utility/libutilityGPU.a
percT.out: /usr/lib/x86_64-linux-gnu/libCGAL.so.13.0.1
percT.out: /usr/local/cuda/lib64/libcudart_static.a
percT.out: /usr/lib/x86_64-linux-gnu/librt.so
percT.out: /usr/lib/x86_64-linux-gnu/libgmp.so
percT.out: /usr/lib/x86_64-linux-gnu/libgmpxx.so
percT.out: /usr/lib/x86_64-linux-gnu/libmpfr.so
percT.out: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
percT.out: /usr/lib/x86_64-linux-gnu/libpthread.so
percT.out: CMakeFiles/percT.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable percT.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/percT.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/percT.out.dir/build: percT.out

.PHONY : CMakeFiles/percT.out.dir/build

CMakeFiles/percT.out.dir/requires: CMakeFiles/percT.out.dir/percT.cpp.o.requires

.PHONY : CMakeFiles/percT.out.dir/requires

CMakeFiles/percT.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/percT.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/percT.out.dir/clean

CMakeFiles/percT.out.dir/depend:
	cd /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build /mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build/CMakeFiles/percT.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/percT.out.dir/depend

