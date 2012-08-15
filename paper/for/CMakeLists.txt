if("${CMAKE_Fortran_COMPILER}" STREQUAL "/usr/bin/gfortran")
  message(FATAL_ERROR "Please specify the Fortran compiler (-DCMAKE_Fortran_COMPILER=gfortran|...)")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -fwhole-file -Ofast")
set(CMAKE_Fortran_FLAGS_DEBUG "-std=f2008 -ffree-form -g -frange-check -Warray-bounds")

add_executable(test-for listing_00.f listing_01.f listing_02.f listing_03.f listing_04.f listing_05.f listing_06.f test.f)
