# load triqs helper to set up tests
set(all_tests
    elk_convert
    elk_equiv_convert
    elk_bands_convert
    elk_bandcharacter_convert
    elk_transport_convert
)

file(GLOB all_test_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
FILE(COPY ${all_test_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${all_tests} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(test_name ${all_tests})
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${test_name}.py)
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()
