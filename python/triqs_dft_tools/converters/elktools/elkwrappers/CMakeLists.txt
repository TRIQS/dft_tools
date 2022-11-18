
# List the sources
set(module_name "getpmatelk")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/getpmatelk.f90") 

set(generated_module_file ${module_name}${TRIQS_PYTHON_MODULE_EXT})

add_custom_target(${module_name} ALL
  DEPENDS ${generated_module_file} 
  )

##generate the fortran python wrapper shared library
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${TRIQS_PYTHON_EXECUTABLE} -m numpy.f2py --f90exec=${CMAKE_Fortran_COMPILER} -c ${fortran_src_file} -m ${module_name} > elk_f2py.log
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

# where to install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} DESTINATION
${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME}/converters/elktools FILES_MATCHING PATTERN "*.so"
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
WORLD_READ WORLD_EXECUTE PATTERN "CMakeFiles" EXCLUDE)
message(STATUS "foo include dir: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "foo include dir: ${CMAKE_BINARY_DIR}")
message(STATUS "foo include dir: ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME}")

# user warning
message(STATUS "-----------------------------------------------------------------------------")
message(STATUS "                             ********  USER NOTE ******** ")
message(STATUS "  This version of DFTTools contains interface routines to read Elk's binary  ")
message(STATUS "  files. ")
message(STATUS "-----------------------------------------------------------------------------")

