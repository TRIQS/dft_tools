# load triqs helper to set up tests
set(TestSuites
    inpconf
#    plocar_io
    plotools
    proj_group
    proj_shell
    vaspio
    atm)
    
FILE(COPY ${TestSuites} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY run_suite.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(test_suite ${TestSuites})
    add_test(${test_suite}
      ${PythonBuildExecutable} run_suite.py ${test_suite})
endforeach(test_suite ${TestSuites})
