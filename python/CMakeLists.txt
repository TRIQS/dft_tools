# where will the python end up in triqs?
set(python_destination applications/dft)

# make a local pytriqs copy
#triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})
triqs_prepare_local_pytriqs(${python_destination})

# to be able to run from toplevel
execute_process(COMMAND ln -fs ${CMAKE_BINARY_DIR}/fortran/vertex/vertex.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})

add_subdirectory(converters/plovasp)
