# where will the python end up in triqs?
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_dft_tools)

# site_customize for build
set(package_name "triqs_dft_tools")

# Create a temporary copy of the python modules so that we can run before installation with the test
FILE(GLOB_RECURSE PYTHON_SOURCES RELATIVE ${CMAKE_SOURCE_DIR}/python *.py )
foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# VASP converter
add_subdirectory(converters/plovasp)

# add version file
configure_file(version.py.in version.py @ONLY)

# install files
install(FILES ${PYTHON_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
