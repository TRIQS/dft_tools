if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT)

 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/app4triqs.modulefile.in ${CMAKE_CURRENT_BINARY_DIR}/app4triqs.modulefile @ONLY)
 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/app4triqs.modulefile DESTINATION ${CMAKE_INSTALL_PREFIX}/share)
 
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/app4triqsvars.sh.in ${CMAKE_CURRENT_BINARY_DIR}/app4triqsvars.sh @ONLY)
 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/app4triqsvars.sh DESTINATION  ${CMAKE_INSTALL_PREFIX}/share)
 
 message(STATUS "***************************************************************")
 message(STATUS "* Custom install Location. Use: 			        ")
 message(STATUS "*                                                              ")
 message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/app4triqsvars.sh      ")
 message(STATUS "*                                                              ") 
 message(STATUS "* to set up the environment variables                          ") 
 message(STATUS "***************************************************************")

endif()
