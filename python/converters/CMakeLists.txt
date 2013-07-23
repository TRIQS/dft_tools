SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_converter.py 
 )

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/dft/converters)
