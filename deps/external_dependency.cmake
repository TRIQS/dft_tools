
function(external_dependency)
  cmake_parse_arguments(ARG "EXCLUDE_FROM_ALL;BUILD_ALWAYS" "VERSION;GIT_REPO;GIT_TAG" "" ${ARGN})

  # -- Was dependency already found?
  get_property(${ARGV0}_FOUND GLOBAL PROPERTY ${ARGV0}_FOUND)
  if(${ARGV0}_FOUND)
    message(STATUS "Dependency ${ARGV0} was already resolved.")
    return()
  endif()

  # -- Try to find package in system.
  if(NOT ARG_BUILD_ALWAYS AND NOT Build_Deps STREQUAL "Always")
    find_package(${ARGV0} ${${ARGV0}_VERSION} QUIET HINTS ${CMAKE_INSTALL_PREFIX})
    if(${ARGV0}_FOUND)
      message(STATUS "Found dependency ${ARGV0} in system.")
      return()
    elseif(Build_Deps STREQUAL "Never")
      message(FATAL_ERROR "Could not find dependency ${ARGV0} in system. Please install the dependency manually or use -DBuild_Deps=IfNotFound during cmake configuration to automatically build all dependencies that are not found.")
    endif()
  endif()

  # -- Build package from source
  message(STATUS " =============== Configuring Dependency ${ARGV0} =============== ")
  if(ARG_EXCLUDE_FROM_ALL)
    set(subdir_opts EXCLUDE_FROM_ALL)
  endif()
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ARGV0})
    message(STATUS "Found sources for dependency ${ARGV0} at ${CMAKE_CURRENT_SOURCE_DIR}/${ARGV0}.")
    add_subdirectory(${ARGV0} ${subdir_opts})
  elseif(ARG_GIT_REPO)
    set(bin_dir ${CMAKE_CURRENT_BINARY_DIR}/${ARGV0})
    set(src_dir ${bin_dir}_src)
    if(NOT IS_DIRECTORY ${src_dir})
      if(ARG_GIT_TAG)
	set(clone_opts --branch ${ARG_GIT_TAG} -c advice.detachedHead=false)
      endif()
      execute_process(COMMAND git clone ${ARG_GIT_REPO} --depth 1 ${clone_opts} ${src_dir})
    endif()
    add_subdirectory(${src_dir} ${bin_dir} ${subdir_opts})
  else()
    message(FATAL_ERROR "Could not find or build dependency ${ARGV0}.")
  endif()

  set_property(GLOBAL PROPERTY ${ARGV0}_FOUND TRUE)

endfunction()
