# List the sources
set(SOURCES modules.f dmftproj.f readcomline.f set_ang_trans.f setsym.f 
            set_rotloc.f timeinv.f read_k_list.f set_projections.f orthogonal.f 
            rot_projectmat.f density.f symmetrize_mat.f rot_dens.f
            orthogonal_wannier.f outputqmc.f outbwin.f outband.f)

# The main target and what to link with...
add_executable(dmftproj ${SOURCES})

# ---------------------------------------------------------------------
# ------ No 1. failed test to get Lapack libs from triqs
# ---------------------------------------------------------------------

# -- trying to TRIQS_LIBRARY_LAPACK from the triqs target
# -- results in ${TRIQS_LIBRARY_LAPACK} == TRIQS_LIBRARY_LAPACK-NOTFOUND
#get_target_property(TRIQS_LIBRARY_LAPACK triqs TRIQS_LIBRARY_LAPACK)

# -- when doing nothing ${TRIQS_LIBRARY_LAPACK} is just empty
message(STATUS "Lapack libraries linked with dmftproj: ${TRIQS_LIBRARY_LAPACK}")

target_link_libraries(dmftproj ${TRIQS_LIBRARY_LAPACK})

# ---------------------------------------------------------------------
# ------ No 2. failed test to get Lapack libs from triqs
# ---------------------------------------------------------------------

#target_link_libraries(dmftproj triqs) # does not work due to "-stdlib=" flag..

# ---------------------------------------------------------------------
# BEGIN Workaround only working on Mac
# ---------------------------------------------------------------------

target_link_libraries(dmftproj "-framework Accelerate")

# ---------------------------------------------------------------------
# END Workaround
# ---------------------------------------------------------------------

# where to install
install (TARGETS dmftproj DESTINATION bin)

# install wien2k files
SET(D ${CMAKE_CURRENT_SOURCE_DIR}/SRC_templates/)
SET(WIEN_SRC_TEMPL_FILES ${D}/case.cf_f_mm2  ${D}/case.cf_p_cubic  ${D}/case.indmftpr ${D}/run_triqs ${D}/runsp_triqs)
message(STATUS "-----------------------------------------------------------------------------")
message(STATUS "                             ********  WARNING ******** ")
message(STATUS "  Wien2k users : after installation of TRIQS, copy the files from ")
message(STATUS "  ${CMAKE_INSTALL_PREFIX}/share/triqs/Wien2k_SRC_files/SRC_templates ")
message(STATUS "  to your Wien2k installation WIENROOT/SRC_templates (Cf documentation).     ")
message(STATUS "  This is not handled automatically by the installation process. ")
message(STATUS "-----------------------------------------------------------------------------")
install (FILES ${WIEN_SRC_TEMPL_FILES} DESTINATION share/triqs/Wien2k_SRC_files/SRC_templates )

