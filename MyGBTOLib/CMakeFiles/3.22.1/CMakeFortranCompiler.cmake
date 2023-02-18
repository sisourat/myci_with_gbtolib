set(CMAKE_Fortran_COMPILER "/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "2021.5.0.20211109")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/intel/oneapi/tbb/2021.5.1/include;/opt/intel/oneapi/mpi/2021.5.1/include;/opt/intel/oneapi/dev-utilities/2021.5.2/include;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/ipp/include;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/include;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/pstl/include;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/tbb/include;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/daal/include;/home/nico/Workspace/Progs/qp2/include;/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/include/intel64;/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/include/icc;/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-linux-gnu/11/include;/usr/include;/usr/include/x86_64-linux-gnu")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/intel/oneapi/tbb/2021.5.1/lib/intel64/gcc4.8;/opt/intel/oneapi/mpi/2021.5.1/libfabric/lib;/opt/intel/oneapi/mpi/2021.5.1/lib/release;/opt/intel/oneapi/mpi/2021.5.1/lib;/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin;/opt/intel/oneapi/compiler/2022.0.2/linux/lib;/opt/intel/oneapi/clck/2021.5.0/lib/intel64;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/mpi/intel64/libfabric/lib;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/ipp/lib/intel64;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/tbb/lib/intel64/gcc4.1;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/daal/lib/intel64_lin;/home/nico/intel/compilers_and_libraries_2019.4.243/linux/tbb/lib/intel64_lin/gcc4.4;/home/nico/Workspace/Progs/qp2/lib;/home/nico/Workspace/Progs/qp2/lib64;/usr/lib/gcc/x86_64-linux-gnu/11;/usr/lib/x86_64-linux-gnu;/usr/lib64;/usr/lib;/lib/x86_64-linux-gnu;/lib64;/lib;/usr/lib/i386-linux-gnu")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
