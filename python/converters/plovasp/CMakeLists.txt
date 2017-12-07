set(python_destination pytriqs/applications/dft/converters/plovasp)

add_cpp2py_module(atm)
target_link_libraries(atm atm_c triqs)
#triqs_set_rpath_for_target(atm)
target_compile_options(atm PRIVATE -std=c++14)

# This we need in order for tests to work
#add_custom_command(TARGET atm POST_BUILD COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}/atm.so ${CMAKE_BINARY_DIR}/python/dft/converters/plovasp)

install (TARGETS atm DESTINATION ${PYTHON_LIB_DEST})


