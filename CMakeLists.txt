# CMake project file for 1WQC_SIM

cmake_minimum_required (VERSION 2.6)
project (1WQC_SIM)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the build directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-llapack -lblas -cpp -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-llapack -lblas -cpp -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-llapack -lblas -cpp -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-llapack -lblas -cpp -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-llapack -lblas -cpp -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-llapack -lblas -cpp -O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#include subdirectories
include_directories ("${PROJECT_SOURCE_DIR}/modules")
add_subdirectory (modules) 

include_directories ("${PROJECT_SOURCE_DIR}/programs")
add_subdirectory (programs) 

# build executables
set (NMPROGRAMS "fivechain" "fivechain_multi" "singlechain" "singlechain_repeating")

foreach (p ${NMPROGRAMS})
  add_executable (${p} "${p}.f")
endforeach (p)

# install executables and scripts
install (TARGETS ${EXECUTABLES} 
         RUNTIME DESTINATION "bin")
