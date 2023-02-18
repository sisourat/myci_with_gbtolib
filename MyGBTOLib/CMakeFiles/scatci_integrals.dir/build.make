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
include CMakeFiles/scatci_integrals.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/scatci_integrals.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/scatci_integrals.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/scatci_integrals.dir/flags.make

CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o: CMakeFiles/scatci_integrals.dir/flags.make
CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o: source/programs/scatci_integrals.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/scatci_integrals.f90 -o CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o

CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.i"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/scatci_integrals.f90 > CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.i

CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.s"
	/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/programs/scatci_integrals.f90 -o CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.s

# Object files for target scatci_integrals
scatci_integrals_OBJECTS = \
"CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o"

# External object files for target scatci_integrals
scatci_integrals_EXTERNAL_OBJECTS =

bin/scatci_integrals: CMakeFiles/scatci_integrals.dir/source/programs/scatci_integrals.f90.o
bin/scatci_integrals: CMakeFiles/scatci_integrals.dir/build.make
bin/scatci_integrals: lib/libGBTO.a
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_core.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/libiomp5.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin/libmkl_core.so
bin/scatci_integrals: /home/nico/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin/libiomp5.so
bin/scatci_integrals: CMakeFiles/scatci_integrals.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable bin/scatci_integrals"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/scatci_integrals.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/scatci_integrals.dir/build: bin/scatci_integrals
.PHONY : CMakeFiles/scatci_integrals.dir/build

CMakeFiles/scatci_integrals.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/scatci_integrals.dir/cmake_clean.cmake
.PHONY : CMakeFiles/scatci_integrals.dir/clean

CMakeFiles/scatci_integrals.dir/depend:
	cd /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/CMakeFiles/scatci_integrals.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/scatci_integrals.dir/depend
