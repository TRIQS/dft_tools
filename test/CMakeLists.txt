# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Copy other files
FILE(COPY SrVO3.pmat SrVO3.struct SrVO3.outputs SrVO3.oubwin SrVO3.ctqmcout SrVO3.symqmc SrVO3.sympar SrVO3.parproj hk_convert_hamiltonian.hk LaVO3-Pnma_hr.dat LaVO3-Pnma.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# List all tests
set(all_tests wien2k_convert hk_convert w90_convert sumkdft_basic srvo3_Gloc srvo3_transp sigma_from_file blockstructure)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath : put the build dir first (in case there is an installed version). 
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )


# VASP converter tests
add_subdirectory(plovasp)

