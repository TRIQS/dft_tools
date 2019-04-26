# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sphinx has internal caching, always run it
add_custom_target(docs_sphinx ALL)
add_custom_command(
  TARGET docs_sphinx
  COMMAND PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

# Autodoc usage requires the python modules to be built first
get_property(CPP2PY_MODULES_LIST GLOBAL PROPERTY CPP2PY_MODULES_LIST)
add_dependencies(docs_sphinx ${CPP2PY_MODULES_LIST})

# -----------------------------------------------------------------------------
# Create an optional target that allows us to regenerate the C++ doc with c++2rst
# -----------------------------------------------------------------------------
foreach(t app4triqs_c triqs)
  get_property(INCLUDE_DIRS TARGET ${t} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  get_property(SYSTEM_INCLUDE_DIRS TARGET ${t} PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
  if(SYTEM_INCLUDE_DIRS)
    list(REMOVE_ITEM INCLUDE_DIRS ${SYSTEM_INCLUDE_DIRS})
  endif()
  foreach(DIR ${INCLUDE_DIRS})
    set(CPP2RST_INCLUDE_FLAGS ${CPP2RST_INCLUDE_FLAGS} -I${DIR})
  endforeach()
  foreach(DIR ${SYSTEM_INCLUDE_DIRS})
    set(CPP2RST_INCLUDE_FLAGS ${CPP2RST_INCLUDE_FLAGS} -isystem=${DIR})
  endforeach()
  get_property(COMPILE_OPTIONS TARGET ${t} PROPERTY INTERFACE_COMPILE_OPTIONS)
  set(CPP2RST_CXXFLAGS ${CPP2RST_CXXFLAGS} ${COMPILE_OPTIONS})
endforeach()

add_custom_target(docs_cpp2rst)
macro(generate_docs header_file)
  add_custom_command(
    TARGET docs_cpp2rst
    COMMAND rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
    COMMAND
      c++2rst
      ${header_file}
      -N app4triqs
      --output_directory ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
      --cxxflags="${CPP2RST_CXXFLAGS}"
      ${CPP2RST_INCLUDE_FLAGS}
  )
endmacro(generate_docs)

generate_docs(${CMAKE_SOURCE_DIR}/c++/app4triqs/toto.hpp)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/app4triqs
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|js|xsl|css|py|txt|inv|bib)$"
  PATTERN "_*"
)
