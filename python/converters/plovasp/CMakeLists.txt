set(python_destination pytriqs/applications/dft/converters/plovasp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})
triqs_python_extension(atm ${python_destination})
target_link_libraries(atm atm_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(atm)

# This we need in order for tests to work
add_custom_command(TARGET atm POST_BUILD COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}/atm.so ${CMAKE_BINARY_DIR}/python/dft/converters/plovasp)

install (TARGETS atm DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
