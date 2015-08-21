# load triqs helper to set up tests
find_package(TriqsTest)

FILE(COPY SrVO3.h5 SrVO3_Sigma.h5 SrVO3.pmat SrVO3.struct SrVO3.outputs SrVO3.oubwin SrVO3.ctqmcout SrVO3.symqmc SrVO3.sympar SrVO3.parproj hk_convert_hamiltonian.hk DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
triqs_add_test_hdf(wien2k_convert " -p 1.e-6" )
triqs_add_test_hdf(hk_convert " -p 1.e-6" )
triqs_add_test_hdf(sumkdft_basic " -d 1.e-6" )
triqs_add_test_hdf(srvo3_Gloc "  -d 1.e-6" )
triqs_add_test_hdf(srvo3_transp " -d 1.e-6" )
triqs_add_test_hdf(sigma_from_file " -d 1.e-6" )
