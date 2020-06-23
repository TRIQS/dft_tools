# List the sources
set(SOURCES modules.f dmftproj.f readcomline.f set_ang_trans.f setsym.f 
            set_rotloc.f timeinv.f read_k_list.f set_projections.f orthogonal.f 
            rot_projectmat.f density.f symmetrize_mat.f rot_dens.f
            orthogonal_wannier.f outputqmc.f outbwin.f outband.f)

# The main target and what to link with...
add_executable(dmftproj ${SOURCES})
target_link_libraries(dmftproj triqs::blas_lapack)

# where to install
install (TARGETS dmftproj DESTINATION bin)

# install wien2k files
SET(D ${CMAKE_CURRENT_SOURCE_DIR}/SRC_templates/)
SET(WIEN_SRC_TEMPL_FILES ${D}/case.cf_f_mm2  ${D}/case.cf_p_cubic  ${D}/case.indmftpr ${D}/run_triqs ${D}/runsp_triqs)
message(STATUS "-----------------------------------------------------------------------------")
message(STATUS "                             ********  WARNING ******** ")
message(STATUS "  Wien2k 14.2 and older : after installation of DFTTools, copy the files from ")
message(STATUS "  ${CMAKE_INSTALL_PREFIX}/share/triqs/Wien2k_SRC_files/SRC_templates ")
message(STATUS "  to your Wien2k installation WIENROOT/SRC_templates (Cf documentation).     ")
message(STATUS "  For newer versions these files are already shipped with Wien2k. ")
message(STATUS "-----------------------------------------------------------------------------")
install (FILES ${WIEN_SRC_TEMPL_FILES} DESTINATION share/triqs/Wien2k_SRC_files/SRC_templates )

