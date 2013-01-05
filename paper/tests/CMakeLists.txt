function(add_tests dir)
  # getting a list of tests
  execute_process(
    COMMAND "find" "${dir}" "-name" "*:in"
    COMMAND "sort"
    COMMAND "cut" "-d" ":" "-f" "1"
    COMMAND "awk" "{printf(\"%s;\",$1)}" # to make it a CMake list
    OUTPUT_VARIABLE TESTS
  )

  # registering the tests to be run by CTest 
  foreach(TEST ${TESTS})
    string(STRIP ${TEST} TEST)
    get_filename_component(NAME "${TEST}" NAME)

    # parsing filenames
    string(REGEX MATCHALL "nx=(.+)_ny=(.+)_Cx=(.+)_Cy=(.+)_nt=(.+)_it=(.+)" NX ${NAME})
    set(NX ${CMAKE_MATCH_1}) # \__ domain size
    set(NY ${CMAKE_MATCH_2}) # /
    set(CX ${CMAKE_MATCH_3}) # \__ Courant numbers
    set(CY ${CMAKE_MATCH_4}) # /
    set(NT ${CMAKE_MATCH_5}) # --- number of timesteps
    set(IT ${CMAKE_MATCH_6}) # --- number of MPDATA iterations
    set(DEC 4) # ----------------- desired precision (abs(desired-actual) < 0.5 * 10**(-decimal))


    set(PFX "--format=time:")
    set(SFX "\;nx=\;${NX}\;nt=\;${NY}\;nt=\;${NT}\;it=\;${IT}\;peakrss=\;%M\;cputime=\;%U\;exit=\;%x\;")


    # C++
    set(TIME_ARG "${PFX}cpp-loac${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
    add_test(NAME "cpp:${NAME}-load" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/cpp/test-cpp ${ARG})
    set_tests_properties("cpp:${NAME}-load" PROPERTIES DEPENDS build_test_cpp)

    set(TIME_ARG "${PFX}cpp-calc${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
    add_test(NAME "cpp:${NAME}-calc" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/cpp/test-cpp ${ARG})
    set_tests_properties("cpp:${NAME}-calc" PROPERTIES DEPENDS build_test_cpp)
  

    # Python
    set(TIME_ARG "${PFX}pyt-load${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
    add_test(NAME "pyt:${NAME}-load" COMMAND time ${TIME_ARG} python ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})

    set(TIME_ARG "${PFX}pyt-calc${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
    add_test(NAME "pyt:${NAME}-calc" COMMAND time ${TIME_ARG} python ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})

    # PyPy
    set(TIME_ARG "${PFX}pyp-load${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
    add_test(NAME "pyp:${NAME}-load" COMMAND time ${TIME_ARG} pypy ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})

    set(TIME_ARG "${PFX}pyp-calc${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
    add_test(NAME "pyp:${NAME}-calc" COMMAND time ${TIME_ARG} pypy ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})

    # Cython
#    set(TIME_ARG "${PFX}cyt${SFX}")
#
#    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
#    add_test(NAME "cyt:${NAME}-load" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/pyt/test-cython ${ARG})
#
#    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
#    add_test(NAME "cyt:${NAME}-calc" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/pyt/test-cython ${ARG})

    # Fortran
    set(TIME_ARG "${PFX}for-load${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
    add_test(NAME "for:${NAME}-load" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/for/test-for ${ARG})
    set_tests_properties("for:${NAME}-load" PROPERTIES DEPENDS build_test_for)

    set(TIME_ARG "${PFX}for-calc${SFX}")
    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
    add_test(NAME "for:${NAME}-calc" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/for/test-for ${ARG})
    set_tests_properties("for:${NAME}-calc" PROPERTIES DEPENDS build_test_for)

    # F77
    #add_test(NAME)
  endforeach()
endfunction()


add_test(build_test_for "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-for)
add_test(build_test_cpp "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-cpp)
add_test(build_test_cython "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-cython)
add_subdirectory(sanity)
add_subdirectory(timing)
