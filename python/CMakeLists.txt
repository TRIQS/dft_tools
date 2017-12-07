# where will the python end up in triqs?
set(python_destination pytriqs/applications/dft)
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/pytriqs/cthyb)

# site_customize for build
set(package_name "pytriqs.applications")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# VASP converter
add_subdirectory(converters/plovasp)

# add version file
configure_file(version.py.in version.py @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})




