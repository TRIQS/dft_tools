
set(python_destination plovasp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})
add_cpp2py_module(atm)
target_link_libraries(atm atm_c)

set(PYTHON_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/atm_desc.py
  ${CMAKE_CURRENT_SOURCE_DIR}/converter.py
  ${CMAKE_CURRENT_SOURCE_DIR}/elstruct.py
  ${CMAKE_CURRENT_SOURCE_DIR}/inpconf.py
  ${CMAKE_CURRENT_SOURCE_DIR}/plotools.py
  ${CMAKE_CURRENT_SOURCE_DIR}/proj_group.py
  ${CMAKE_CURRENT_SOURCE_DIR}/proj_shell.py
  ${CMAKE_CURRENT_SOURCE_DIR}/sc_dmft.py
  ${CMAKE_CURRENT_SOURCE_DIR}/vaspio.py)

set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_dft_tools)
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST}/converters/plovasp)
install (TARGETS atm DESTINATION ${PYTHON_LIB_DEST}/converters/plovasp)
