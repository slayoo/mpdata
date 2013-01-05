function(add_test_helper)
  #TODO
endfunction()

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

    set(REPEAT "0;1;2;")

    # C++
    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}cpp-load-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
      add_test(NAME "cpp:${NAME}-load-${N}" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/cpp/test-cpp ${ARG})
      set_tests_properties("cpp:${NAME}-load-${N}" PROPERTIES DEPENDS build_test_cpp)
    endforeach()

    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}cpp-calc-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
      add_test(NAME "cpp:${NAME}-calc-${N}" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/cpp/test-cpp ${ARG})
      set_tests_properties("cpp:${NAME}-calc-${N}" PROPERTIES DEPENDS build_test_cpp)
    endforeach()
  

    # Python
    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}pyt-load-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
      add_test(NAME "pyt:${NAME}-load-${N}" COMMAND time ${TIME_ARG} python ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})
    endforeach()

    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}pyt-calc-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
      add_test(NAME "pyt:${NAME}-calc-${N}" COMMAND time ${TIME_ARG} python ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})
    endforeach()

    # PyPy
    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}pyp-load-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
      add_test(NAME "pyp:${NAME}-load-${N}" COMMAND time ${TIME_ARG} pypy ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})
    endforeach()

    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}pyp-calc-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
      add_test(NAME "pyp:${NAME}-calc-${N}" COMMAND time ${TIME_ARG} pypy ${CMAKE_BINARY_DIR}/pyt/test.py ${ARG})
    endforeach()

    # Cython
#    set(TIME_ARG "${PFX}cyt${SFX}")
#
#    set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
#    add_test(NAME "cyt:${NAME}-load" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/pyt/test-cython ${ARG})
#
#    set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
#    add_test(NAME "cyt:${NAME}-calc" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/pyt/test-cython ${ARG})

    # Fortran
    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}for-load-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
      add_test(NAME "for:${NAME}-load-${N}" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/for/test-for ${ARG})
      set_tests_properties("for:${NAME}-load-${N}" PROPERTIES DEPENDS build_test_for)
    endforeach()

    foreach(N ${REPEAT})
      set(TIME_ARG "${PFX}for-calc-${N}${SFX}")
      set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
      add_test(NAME "for:${NAME}-calc-${N}" COMMAND time ${TIME_ARG} ${CMAKE_BINARY_DIR}/for/test-for ${ARG})
      set_tests_properties("for:${NAME}-calc-${N}" PROPERTIES DEPENDS build_test_for)
    endforeach()

    # F77
    #add_test(NAME)
  endforeach()
endfunction()


add_test(build_test_for "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-for)
add_test(build_test_cpp "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-cpp)
add_test(build_test_cython "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test-cython)
add_subdirectory(sanity)
add_subdirectory(timing)
