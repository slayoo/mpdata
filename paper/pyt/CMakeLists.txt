#find_package(PythonLibs)
#include_directories(${PYTHON_INCLUDE_DIRS})

#add_custom_target(test-c cython -X boundscheck=False --embed test.py)
#add_dependencies(test-c test.py)

#add_custom_target(listings-c cython -X boundscheck=False listings.py)
#add_dependencies(listings-c listings.py)

#add_executable(test-cython test.c listings.c)
#add_dependencies(test-cython test-c listings-c)
#target_link_libraries(test-cython ${PYTHON_LIBRARIES})
