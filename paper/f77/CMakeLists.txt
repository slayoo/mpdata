set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffree-form -fwhole-file -Ofast")
set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -ffree-form -g -frange-check -Warray-bounds")

add_executable(test-f77-3x3 test.f)
set_target_properties(test-f77-3x3 PROPERTIES COMPILE_FLAGS "-DNX=3 -DNY=3")

add_executable(test-f77-64x64 test.f)
set_target_properties(test-f77-64x64 PROPERTIES COMPILE_FLAGS "-DNX=64 -DNY=64")
