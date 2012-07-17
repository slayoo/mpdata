cmake_minimum_required(VERSION 2.6) # for simplified if-constructs

if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "Release"))
  message(FATAL_ERROR "Please specify the build type, i.e.:
    \$ cmake . -DCMAKE_BUILD_TYPE=Debug
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=Release 
  ")
endif()

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "Makefile;CMakeCache.txt;cmake_install.cmake;CMakeFiles;")

project(egu2012 CXX Fortran)

### Fortran stuff:
add_executable(egu2012-for egu2012.f)
target_link_libraries(egu2012-for "-lnetcdff")

set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -cpp -Ofast -fwhole-file" CACHE STRING "Release options." FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-std=f2008 -ffree-form -cpp -g -Warray-bounds" CACHE STRING "Debug options." FORCE)

### C++ stuff:
add_executable(egu2012-cpp egu2012.cpp)

#set(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -Ofast -DNDEBUG" CACHE STRING "Release options." FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -DNDEBUG" CACHE STRING "Release options." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x -g -DBZ_DEBUG" CACHE STRING "Debug options." FORCE)

include(CheckCXXSourceCompiles) 
include(CheckCXXSourceRuns)

# Blitz
check_cxx_source_compiles( 
  "#include <blitz/array.h>\n int main() {}"
  Find_Blitz
)
if(NOT Find_Blitz) 
  message(FATAL_ERROR "Blitz++ not found!
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()                                                                
check_cxx_source_runs("
  #include <blitz/array.h>
  int main() { 
    assert(std::string(BZ__compiler_name) == std::string(\"${CMAKE_CXX_COMPILER}\"));
  }
  " Blitz_Configured_With_The_Same_Compiler
)
if (NOT Blitz_Configured_With_The_Same_Compiler)
  message(WARNING "Blitz++ seem not to have been configured with the same compiler!
    Consider reconfiguring Blitz++ with CXX=${CMAKE_CXX_COMPILER} 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif() 
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
  target_link_libraries(egu2012-cpp "-lblitz")
endif()

# OpenMP (Not using CMake's FIND_PACKAGE(OpenMP) as it relies on C compiler rather than C++ one)
if ("${OPENMP}" STREQUAL "NO")   
  message(STATUS "Building without OpenMP (-DOPENMP=NO specified)")
else()
  set(TMP ${CMAKE_CXX_FLAGS}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  check_cxx_source_compiles(
    "#include <omp.h>\n int main() { omp_set_num_threads(2); }" 
    Find_OpenMP
  )
  if(NOT Find_OpenMP)
    message(WARNING "OpenMP not found (cf. CMakeFiles/CMakeError.log) 
      OpenMP parallelisation will not be available.")
    set(CMAKE_CXX_FLAGS ${TMP}) 
  endif()   
endif()

#NetCDF
set(CMAKE_REQUIRED_LIBRARIES "-lnetcdf_c++4")
check_cxx_source_compiles(
  "#include <netcdf>\n int main() { netCDF::NcFile *f; }" 
  Find_NetCDF
)
if(NOT Find_NetCDF) 
  message(FATAL_ERROR "
      NetCDF >=4.1.1 not found! (cf. CMakeFiles/CMakeError.log)
      Without it the NetCDF output will not be available.
      After installing it erase the CMakeCache.txt before re-running cmake
    ")
else()
  set(USE_NETCDF ${Find_NetCDF})
  target_link_libraries(egu2012-cpp ${CMAKE_REQUIRED_LIBRARIES})
endif() 
