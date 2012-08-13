cmake_minimum_required(VERSION 2.6.3)
project(mpdata CXX Fortran)

if(NOT (
  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" 
))
  message(FATAL_ERROR "Please specify the build type (-DCMAKE_BUILD_TYPE=Debug|Release") 
endif()

add_subdirectory(cpp)
add_subdirectory(for)
