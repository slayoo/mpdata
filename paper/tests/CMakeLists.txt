execute_process(
  COMMAND "find" "${CMAKE_SOURCE_DIR}/tests" "-name" "*:in"
  COMMAND "cut" "-d" ":" "-f" "1"
  OUTPUT_VARIABLE TESTS
)

foreach(TEST ${TESTS})
  string(STRIP ${TEST} TEST)

  get_filename_component(NAME "${TEST}" NAME)
  string(REGEX MATCHALL "nx=(.+)_ny=(.+)_Cx=(.+)_Cy=(.+)_nt=(.+)_it=(.+)" NX ${NAME})
  set(NX ${CMAKE_MATCH_1})
  set(NY ${CMAKE_MATCH_2})
  set(CX ${CMAKE_MATCH_3})
  set(CY ${CMAKE_MATCH_4})
  set(NT ${CMAKE_MATCH_5})
  set(IT ${CMAKE_MATCH_6})
  set(DEC 5)

  # C++
  add_test(NAME "cpp:${NAME}" COMMAND ${CMAKE_BINARY_DIR}/cpp/test-cpp 
    ${NX} ${NY} ${CX} ${CY} ${NT} ${IT} ${NAME}:in ${NAME}:out ${DEC}
  )
  
  # Python
  add_test(NAME "pyt:${NAME}" COMMAND python ${CMAKE_BINARY_DIR}/pyt/test.py 
    ${NX} ${NY} ${CX} ${CY} ${NT} ${IT} ${NAME}:in ${NAME}:out ${DEC}
  )

  # TODO: PyPy

  # Fortran
  add_test(NAME "for:${NAME}" COMMAND ${CMAKE_BINARY_DIR}/for/test-for 
    ${NX} ${NY} ${CX} ${CY} ${NT} ${IT} ${NAME}:in ${NAME}:out ${DEC}
  )

endforeach()
