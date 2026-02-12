# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
if(NOT DEFINED SPHINXBUILD_EXECUTABLE)
  find_package(Sphinx)
endif()

# Sphinx has internal caching, always run it
add_custom_target(${PROJECT_NAME}_docs_sphinx ALL)
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx
  POST_BUILD
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:${triqs_dftkit_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SPHINXBUILD_EXECUTABLE} -j auto -c . -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/${PROJECT_NAME}
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|ico|js|xsl|css|py|txt|inv|bib|ttf|woff2|eot|INCAR|POSCAR|cfg|KPOINTS)$"
  PATTERN "_*"
)
