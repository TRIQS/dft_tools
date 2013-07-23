# List the sources
set (SOURCES modules dmftproj readcomline set_ang_trans setsym 
     set_rotloc timeinv read_k_list set_projections orthogonal 
     rot_projectmat density symmetrize_mat rot_dens
     orthogonal_wannier outputqmc outbwin outband)

# add the extension and the path
FOREACH(f ${SOURCES} )
 set(S "${CMAKE_CURRENT_SOURCE_DIR}/${f}.f;${S}")
ENDFOREACH(f)
   
# The main target and what to link with...
add_executable(dmftproj ${S})
target_link_libraries(dmftproj ${TRIQS_LIBRARY_LAPACK} )

# where to install
install (TARGETS dmftproj DESTINATION bin )

# that is it !

SET( D ${CMAKE_CURRENT_SOURCE_DIR}/SRC_templates/)
SET(WIEN_SRC_TEMPL_FILES ${D}/case.cf_f_mm2  ${D}/case.cf_p_cubic  ${D}/case.indmftpr ${D}/run_triqs ${D}/runsp_triqs)

# build the fortran stuff...
message(STATUS "-----------------------------------------------------------------------------")
message(STATUS "                             ********  WARNING ******** ")
message(STATUS "  Wien2k users : after installation of TRIQS, copy the files from ")
message(STATUS "  ${CMAKE_INSTALL_PREFIX}/share/triqs/Wien2k_SRC_files/SRC_templates ")
message(STATUS "  to your Wien2k installation WIENROOT/SRC_templates (Cf documentation).     ")
message(STATUS "  This is not handled automatically by the installation process. ")
message(STATUS "-----------------------------------------------------------------------------")
install (FILES ${WIEN_SRC_TEMPL_FILES} DESTINATION share/triqs/Wien2k_SRC_files/SRC_templates )

