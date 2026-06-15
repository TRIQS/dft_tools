# Fortran regression tests for the dmftproj symmetry construction.
# Compile the dmftproj sources (excluding the main program) together with the
# driver so the real dmat / d_matrix / timeinv_op are exercised. See #148.
file(GLOB dmftproj_sources ${PROJECT_SOURCE_DIR}/fortran/dmftproj/*.f)
list(REMOVE_ITEM dmftproj_sources ${PROJECT_SOURCE_DIR}/fortran/dmftproj/dmftproj.f)

add_executable(test_spinrot_symmetry test_spinrot_symmetry.f90 ${dmftproj_sources})
target_link_libraries(test_spinrot_symmetry nda::blas_lapack)
add_test(NAME test_spinrot_symmetry COMMAND test_spinrot_symmetry)
