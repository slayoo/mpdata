find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})

add_custom_target(test-c cython --embed test.py)
add_custom_target(listings-c cython listings.py)
add_executable(test-cython test.c listings.c)
add_dependencies(test-cython test-c listings-c)
target_link_libraries(test-cython ${PYTHON_LIBRARIES})
