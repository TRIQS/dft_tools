
FILE(GLOB TestList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 target_link_libraries( ${TestName} atm_c ${TRIQS_LIBRARY_ALL} )

 add_test(NAME ${TestName} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TestName})

 if (TESTS_C_WITH_VALGRIND)
  add_test ( ${TestName}_valgrind valgrind --error-exitcode=1 ${CMAKE_CURRENT_BINARY_DIR}/${TestName})
 endif()

ENDFOREACH( TestName1  ${TestList} )
