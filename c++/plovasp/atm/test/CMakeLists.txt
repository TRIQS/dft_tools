enable_testing()

FILE(GLOB TestList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}  ${TestName}.cpp )
 target_link_libraries( ${TestName} atm_c triqs)
 triqs_set_rpath_for_target( ${TestName} )
 add_test(NAME ${TestName} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${t})
ENDFOREACH( TestName1  ${TestList} )



