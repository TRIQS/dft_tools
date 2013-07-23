
triqs_make_target_to_copy_all_py_files_from_python_dir_to_build_dir()

execute_process(COMMAND ln -fs ${CMAKE_BINARY_DIR}/fortran/F90/vertex.so ${CMAKE_CURRENT_BINARY_DIR} )

add_subdirectory(converters)

#installation 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/solver_multiband.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/sumk_lda.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/sumk_lda_tools.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/U_matrix.py 
 )
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/dft)

