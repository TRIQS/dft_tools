# load triqs helper to set up tests
find_package(TriqsTest)

FILE(COPY SrVO3.h5 SrVO3.ctqmcout SrVO3.symqmc SrVO3.sympar SrVO3.parproj hamiltonian.hk DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
triqs_add_test_hdf(wien2k_convert " -p 1.e-6" )
triqs_add_test_hdf(hk_convert " -p 1.e-6" )
triqs_add_test_hdf(sumklda_basic " -d 1.e-6" )
triqs_add_test_hdf(srvo3_Gloc "  -d 1.e-6" )
triqs_add_test_hdf(U_mat " -d 1.e-6" )

