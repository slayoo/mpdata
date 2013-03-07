cmake_minimum_required(VERSION 2.8.0) 
project(mpdata C CXX Fortran)

#if(NOT (
#  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
#  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" 
#))
#  message(FATAL_ERROR "Please specify the build type (-DCMAKE_BUILD_TYPE=Debug|Release)") 
#endif()

enable_testing()
add_subdirectory(cpp)
add_subdirectory(for)
add_subdirectory(pyt)
#add_subdirectory(f77)
add_subdirectory(tests)
