# where will the python end up in triqs?
set(python_destination pytriqs/applications/dft)
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/pytriqs/cthyb)

# site_customize for build
set(package_name "pytriqs.applications")

# Create a temporary copy of the python modules so that we can run before installation with the test
FILE(GLOB_RECURSE all_py_files RELATIVE ${CMAKE_SOURCE_DIR}/python *.py )
foreach(f ${all_py_files})
 configure_file(${f} ${f} COPYONLY)
endforeach()

#configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# VASP converter
add_subdirectory(converters/plovasp)

# add version file
configure_file(version.py.in version.py @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})




