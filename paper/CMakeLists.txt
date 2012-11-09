cmake_minimum_required(VERSION 2.6.3) # TODO: "module procedure bug" http://www.cmake.org/pipermail/cmake/2007-October/016831.html
project(mpdata CXX Fortran)

if(NOT (
  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" 
))
  message(FATAL_ERROR "Please specify the build type (-DCMAKE_BUILD_TYPE=Debug|Release)") 
endif()

enable_testing()
add_subdirectory(cpp)
add_subdirectory(for)
add_subdirectory(tests)
