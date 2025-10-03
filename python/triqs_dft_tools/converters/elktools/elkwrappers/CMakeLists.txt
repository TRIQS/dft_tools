# F2PY headers
execute_process(
  COMMAND "${TRIQS_PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT getpmatelkmodule.c getpmatelk-f2pywrappers.f getpmatelk-f2pywrappers2.f90
  DEPENDS getpmatelk.f90
  VERBATIM
  COMMAND "${CMAKE_COMMAND}" -E touch "getpmatelkmodule.c" "getpmatelk-f2pywrappers.f" "getpmatelk-f2pywrappers2.f90"
  COMMAND "${TRIQS_PYTHON_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/getpmatelk.f90" -m getpmatelk --lower)

Python_add_library(getpmatelk MODULE
                   "${CMAKE_CURRENT_BINARY_DIR}/getpmatelkmodule.c"
                   "${CMAKE_CURRENT_BINARY_DIR}/getpmatelk-f2pywrappers.f"
                   "${CMAKE_CURRENT_BINARY_DIR}/getpmatelk-f2pywrappers2.f90"
                   "${CMAKE_CURRENT_SOURCE_DIR}/getpmatelk.f90")
set_property(TARGET getpmatelk PROPERTY SUFFIX "${TRIQS_PYTHON_MODULE_EXT}")
target_link_libraries(getpmatelk PRIVATE fortranobject)

install(TARGETS getpmatelk DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME}/converters/elktools/elkwrappers)

# user warning
message(STATUS "-----------------------------------------------------------------------------")
message(STATUS "                             ********  USER NOTE ******** ")
message(STATUS "  This version of DFTTools contains interface routines to read Elk's binary  ")
message(STATUS "  files. ")
message(STATUS "-----------------------------------------------------------------------------")
