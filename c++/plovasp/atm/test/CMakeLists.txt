find_package(TriqsTest)
enable_testing()

# Linking and include info
#add_library(atm_c dos_tetra3d.hpp dos_tetra3d.cpp argsort.h argsort.c)
#set_target_properties(atm_c PROPERTIES LINKER_LANGUAGE CXX)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/c++/plovasp/atm ${TRIQS_INCLUDE_ALL})

FILE(GLOB TestList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 target_link_libraries( ${TestName} atm_c ${TRIQS_LIBRARY_ALL} )
 triqs_set_rpath_for_target( ${TestName} )
 triqs_add_cpp_test( ${TestName}    )
# if (TESTS_C_WITH_VALGRIND)
#  add_test ( ${TestName}_valgrind valgrind --error-exitcode=1 ${CMAKE_CURRENT_BINARY_DIR}/${TestName})
# endif()
ENDFOREACH( TestName1  ${TestList} )

#add_executable(test_atm test2py.cpp)
#target_link_libraries(test_atm atm_c)

#add_subdirectory(test)
