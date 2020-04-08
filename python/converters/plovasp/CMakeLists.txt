# === Build and install atm module 
add_cpp2py_module(atm)
target_link_libraries(atm PRIVATE atm_c triqs)
target_compile_options(atm PRIVATE -std=c++17)
target_include_directories(atm PRIVATE ${CMAKE_SOURCE_DIR}/c++)

install(TARGETS atm DESTINATION ${PYTHON_LIB_DEST}/converters/plovasp)

# === Copy Python files to current build directory and register for install
set(PYTHON_SOURCES __init__.py  converter.py elstruct.py  inpconf.py  plotools.py  proj_group.py  proj_shell.py  sc_dmft.py  vaspio.py)
foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# install files
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST}/converters/plovasp)

# This we need in order for tests to work
#add_custom_command(TARGET atm POST_BUILD COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}/atm.so ${CMAKE_BINARY_DIR}/python/dft/converters/plovasp)
