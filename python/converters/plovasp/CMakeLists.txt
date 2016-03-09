include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})
triqs_python_extension(atm)
target_link_libraries(atm atm_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(atm)

