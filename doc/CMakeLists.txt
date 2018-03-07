# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

set(EXT_FOR_DOC "rst png txt css_t conf css js gif jpg py html bib sh")
execute_process(COMMAND cp_rs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${EXT_FOR_DOC})

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sources
file(GLOB_RECURSE sources *.rst)

# create documentation target
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${sources}
  		   COMMAND ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_BINARY_DIR} html)
add_custom_target(doc_sphinx ALL DEPENDS ${sphinx_top} ${CMAKE_CURRENT_BINARY_DIR})

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/dft_tools
 FILES_MATCHING
 REGEX "\\.(html|pdf|png|gif|jpg|js|xsl|css|py|txt|inv|bib)$"
 PATTERN "_*"
 )
