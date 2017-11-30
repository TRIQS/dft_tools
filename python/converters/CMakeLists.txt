
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_dft_tools)

set(PYTHON_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/converter_tools.py
  ${CMAKE_CURRENT_SOURCE_DIR}/hk_converter.py
  ${CMAKE_CURRENT_SOURCE_DIR}/vasp_converter.py
  ${CMAKE_CURRENT_SOURCE_DIR}/wannier90_converter.py
  ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_converter.py)

install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST}/converters)

# VASP converter
add_subdirectory(plovasp)
