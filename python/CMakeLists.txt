# site_customize for build
#set(package_name "pytriqs.applications")

#configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# The std for all targets
add_compile_options(-std=c++1z)

# add version file
configure_file(version.py.in version.py)

set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_dft_tools)

set(PYTHON_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/block_structure.py
  ${CMAKE_CURRENT_SOURCE_DIR}/clear_h5_output.py
  ${CMAKE_CURRENT_SOURCE_DIR}/sumk_dft.py
  ${CMAKE_CURRENT_SOURCE_DIR}/sumk_dft_tools.py
  ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.py
  ${CMAKE_CURRENT_SOURCE_DIR}/trans_basis.py
  ${CMAKE_CURRENT_SOURCE_DIR}/update_archive.py)

install(FILES ${PYTHON_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})

add_subdirectory(converters)
