
# List the sources
set(module_name "getpmatelk")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/getpmatelk.f90") 

execute_process(COMMAND ${TRIQS_PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or sysconfig.get_config_var('SO'))"
              OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
              RESULT_VARIABLE FOUND_PYTHON_EXT_SUFFIX
              OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)
message(FATAL_ERROR "Unable to determine file extension of compiled Python modules - specify it with F2PY_SUFFIX")
endif (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)

message(STATUS ${module_name}${PYTHON_EXT_SUFFIX})

set(generated_module_file ${module_name}${PYTHON_EXT_SUFFIX})

add_custom_target(${module_name} ALL
  DEPENDS ${generated_module_file} 
  )


##generate the fortran python wrapper shared library
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${TRIQS_PYTHON_EXECUTABLE} -m numpy.f2py --f90exec=${CMAKE_Fortran_COMPILER} --noarch -c ${fortran_src_file} -m ${module_name} > elk_f2py.log
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

