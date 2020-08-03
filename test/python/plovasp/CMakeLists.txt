# load triqs helper to set up tests
set(all_tests
    atm
    converter
    inpconf
    plotools
    proj_group
    proj_shell
    vaspio
)
    
FILE(COPY ${all_tests} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY run_suite.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(test_name ${all_tests})
  add_test(NAME Py_vasp_${test_name} COMMAND ${TRIQS_PYTHON_EXECUTABLE} run_suite.py ${test_name})
  set_property(TEST Py_vasp_${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()


