if("${CMAKE_Fortran_COMPILER}" STREQUAL "/usr/bin/gfortran")
  message(FATAL_ERROR "Please specify the Fortran compiler (-DCMAKE_Fortran_COMPILER=gfortran|...)")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -cpp -fwhole-file -Ofast")
set(CMAKE_Fortran_FLAGS_DEBUG "-std=f2008 -ffree-form -cpp -g -frange-check -Warray-bounds")

add_executable(test-for test.f)
