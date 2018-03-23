# Generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sources
file(GLOB_RECURSE sources *.rst)

set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/contents.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${sources}
  		   COMMAND ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_BINARY_DIR} html)
add_custom_target(docs_sphinx ALL DEPENDS ${sphinx_top} ${CMAKE_CURRENT_BINARY_DIR})

# Application must be compiled first
add_dependencies(docs_sphinx app4triqs_c)

# ------------------------------------------------------------------------------------------------
# Copy the directory substructure and link all files with relevant extensions
# cp_rs is a script in cpp2py/bin, it mimics cp -rs on Linux
# and filters the relevant extension
# ------------------------------------------------------------------------------------------------

set(EXT_FOR_DOC "rst png txt css_t conf css js gif jpg py html bib")
execute_process(COMMAND cp_rs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${EXT_FOR_DOC})

# ---------------------------------
# Generate C++ doc with c++2rst
# ---------------------------------

add_custom_target(docs_cpp2rst)

get_property(TRIQS_INCLUDE_DIRS TARGET triqs PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(APP4TRIQS_INCLUDE_DIRS TARGET app4triqs_c PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

set(CPP2RST_INCLUDE_DIRS ${TRIQS_INCLUDE_DIRS} ${APP4TRIQS_INCLUDE_DIRS})
foreach(I ${CPP2RST_INCLUDE_DIRS})
  set (CPP2RST_INCLUDE_COMMAND ${CPP2RST_INCLUDE_COMMAND} --includes ${I})
endforeach()

macro(generate_docs header_file)
  add_custom_command (TARGET docs_cpp2rst
			COMMAND c++2rst
			${header_file}
			-N app4triqs
			--output_directory ${CMAKE_CURRENT_BINARY_DIR}/cpp2rst_generated
			${CPP2RST_INCLUDE_COMMAND}
			--cxxflags="-std=c++17"
			-I ${CMAKE_SOURCE_DIR}/c++
			-I ${TRIQS_ROOT}/include
			2>&1 >> cpp2rst.log
		      )
endmacro(generate_docs)

generate_docs(${CMAKE_SOURCE_DIR}/c++/app4triqs/toto.hpp)

add_dependencies(docs_sphinx docs_cpp2rst)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/app4triqs
 FILES_MATCHING
 REGEX "\\.(html|pdf|png|gif|jpg|js|xsl|css|py|txt|inv|bib)$"
 PATTERN "_*"
 )

