triqs_build_f2py_module( triqs_DFT vertex vertex.pyf vertex.f90)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/vertex.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/dft)
