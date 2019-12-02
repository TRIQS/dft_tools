# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# -----------------------------------------------------------------------------
# Create an optional target that allows us to regenerate the C++ doc with c++2rst
# -----------------------------------------------------------------------------
add_custom_target(docs_cpp2rst)
include(${PROJECT_SOURCE_DIR}/share/cmake/extract_flags.cmake)
extract_flags(app4triqs_c)
separate_arguments(app4triqs_c_CXXFLAGS)
macro(generate_docs header_file)
  add_custom_command(
    TARGET docs_cpp2rst
    COMMAND rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
    COMMAND
      c++2rst
      ${header_file}
      -N app4triqs
      --output_directory ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
      -I${PROJECT_SOURCE_DIR}/c++
      --cxxflags="${app4triqs_c_CXXFLAGS}"
  )
endmacro(generate_docs)

generate_docs(${PROJECT_SOURCE_DIR}/c++/app4triqs/app4triqs.hpp)

# --------------------------------------------------------
# Build & Run the C++ doc examples and capture the output
# --------------------------------------------------------

add_custom_target(docs_example_output)
file(GLOB_RECURSE ExampleList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach(example ${ExampleList})
  get_filename_component(f ${example} NAME_WE)
  get_filename_component(d ${example} DIRECTORY)
  add_executable(doc_${f} EXCLUDE_FROM_ALL ${example})
  set_property(TARGET doc_${f} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${d})
  target_link_libraries(doc_${f} triqs)
  add_custom_command(TARGET doc_${f}
   COMMAND doc_${f} > ${CMAKE_CURRENT_SOURCE_DIR}/${d}/${f}.output 2>/dev/null
   WORKING_DIRECTORY ${d}
  )
  add_dependencies(docs_example_output doc_${f})
endforeach()

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sphinx has internal caching, always run it
add_custom_target(docs_sphinx ALL)
add_custom_command(
  TARGET docs_sphinx
  COMMAND PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

option(Sphinx_Only "When building the documentation, skip the Python Modules and the generation of C++ Api and example outputs" OFF)
if(NOT Sphinx_Only)
  # Autodoc usage requires the python modules to be built first
  get_property(CPP2PY_MODULES_LIST GLOBAL PROPERTY CPP2PY_MODULES_LIST)
  if(CPP2PY_MODULES_LIST)
    add_dependencies(docs_sphinx ${CPP2PY_MODULES_LIST})
  endif()

  # Generation of C++ Api and Example Outputs
  add_dependencies(docs_sphinx docs_cpp2rst docs_example_output)
endif()

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/app4triqs
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib)$"
  PATTERN "_*"
)
