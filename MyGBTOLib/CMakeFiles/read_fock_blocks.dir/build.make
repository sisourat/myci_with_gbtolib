# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib

# Include any dependencies generated for this target.
include CMakeFiles/read_fock_blocks.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/read_fock_blocks.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/read_fock_blocks.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/read_fock_blocks.dir/flags.make

CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o: CMakeFiles/read_fock_blocks.dir/flags.make
CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o: source/programs/read_fock_blocks.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/read_fock_blocks.f90 -o CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o

CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.i"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/read_fock_blocks.f90 > CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.i

CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.s"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/read_fock_blocks.f90 -o CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.s

# Object files for target read_fock_blocks
read_fock_blocks_OBJECTS = \
"CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o"

# External object files for target read_fock_blocks
read_fock_blocks_EXTERNAL_OBJECTS =

bin/read_fock_blocks: CMakeFiles/read_fock_blocks.dir/source/programs/read_fock_blocks.f90.o
bin/read_fock_blocks: CMakeFiles/read_fock_blocks.dir/build.make
bin/read_fock_blocks: lib/libGBTO.a
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_core.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/libiomp5.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_core.so
bin/read_fock_blocks: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/libiomp5.so
bin/read_fock_blocks: CMakeFiles/read_fock_blocks.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable bin/read_fock_blocks"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/read_fock_blocks.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/read_fock_blocks.dir/build: bin/read_fock_blocks
.PHONY : CMakeFiles/read_fock_blocks.dir/build

CMakeFiles/read_fock_blocks.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/read_fock_blocks.dir/cmake_clean.cmake
.PHONY : CMakeFiles/read_fock_blocks.dir/clean

CMakeFiles/read_fock_blocks.dir/depend:
	cd /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles/read_fock_blocks.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/read_fock_blocks.dir/depend

