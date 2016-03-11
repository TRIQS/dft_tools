set(python_destination pytriqs/applications/dft/converters/plovasp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})
triqs_python_extension(atm ${python_destination})
target_link_libraries(atm atm_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(atm)

#execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/atm.so ${CMAKE_BINARY_DIR}/${python_destination})

install (TARGETS atm DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
