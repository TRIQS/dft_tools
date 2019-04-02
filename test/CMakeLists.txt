add_subdirectory(c++)

if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
endif()
