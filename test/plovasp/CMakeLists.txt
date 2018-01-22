# load triqs helper to set up tests
set(all_tests
    inpconf
#    plocar_io
    plotools
    proj_group
    proj_shell
    vaspio
    atm)
    
FILE(COPY ${all_tests} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY run_suite.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python run_suite.py ${t})
endforeach()

set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
