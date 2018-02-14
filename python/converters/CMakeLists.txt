# Create a temporary copy of the python modules so that we can run before installation with the test
FILE(GLOB PYTHON_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# install files
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST}/converters)

add_subdirectory(plovasp)
