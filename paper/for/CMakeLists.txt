set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -fwhole-file -Ofast -ftree-vectorizer-verbose=1")
set(CMAKE_Fortran_FLAGS_DEBUG "-std=f2008 -fcheck=all -ffree-form -g -Warray-bounds -ftree-vectorizer-verbose=1")

add_executable(test-for listings.f test.f)
add_executable(test_arrvec listings.f test_arrvec.f)
add_test(test_arrvec test_arrvec)
