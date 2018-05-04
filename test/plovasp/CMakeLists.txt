# load triqs helper to set up tests
set(all_tests
    inpconf
    proj_group
    proj_shell
    vaspio
    atm
    plotools
    converter
)
    
FILE(COPY ${all_tests} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY run_suite.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND ${python_executable} run_suite.py ${t})
endforeach()

set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
