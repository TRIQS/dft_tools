# Build the python module
add_cpp2py_module(toto_module)
target_link_libraries(toto_module app4triqs_c)

# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
list(FILTER python_sources EXCLUDE REGEX "_desc.py")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/app4triqs)
install(TARGETS toto_module DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${python_sources} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
