# where will the python end up in triqs?
set(python_destination pytriqs/applications/dft)

# site_customize for build
set(package_name "pytriqs.applications")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# make a local pytriqs copy
triqs_prepare_local_pytriqs(${python_destination})
