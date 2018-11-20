# Start configuration
cmake_minimum_required(VERSION 2.8.12)
project(app4triqs CXX)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Build static libraries
set(BUILD_SHARED_LIBS OFF)

# Enable compiler warnings for the whole project
add_compile_options(-Wall
		     $<$<CONFIG:Debug>:-Og>
		     $<$<CONFIG:Debug>:-ggdb3>
		     )

# Load Dependencies
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Use git to detect app4triqs version and git hash
execute_process(COMMAND bash "-c" "${TRIQS_GIT_EXECUTABLE} describe | sed \"s/-\\([0-9]*\\)-.*/.\\1/g\""
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   APP4TRIQS_VERSION
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT APP4TRIQS_VERSION)
  set(APP4TRIQS_VERSION 0.1)
endif()
triqs_get_git_hash_of_source_dir(APP4TRIQS_GIT_HASH)
MESSAGE(STATUS "app4triqs version : ${APP4TRIQS_VERSION}")
message(STATUS "Git hash: ${APP4TRIQS_GIT_HASH}")

# Build and install the app4triqs library
add_subdirectory(c++/app4triqs)

# Build and install the app4triqs python module
if(${TRIQS_WITH_PYTHON_SUPPORT})
  add_subdirectory(python/app4triqs)
endif()

# Tests
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
  add_subdirectory(test)
endif()

# Additional configuration files
add_subdirectory(share)

option(Build_Documentation "Build documentation" OFF)
if(Build_Documentation)
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message(WARNING "TRIQS library has not been compiled with its documentation! Cannot build documentation.")
  else()
    message(STATUS "Documentation Build enabled")
    add_subdirectory(doc)
  endif()
endif()

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  SET(CPACK_GENERATOR "DEB")
  SET(CPACK_PACKAGE_VERSION ${APP4TRIQS_VERSION})
  SET(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/app4triqs")
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 1.5), cpp2py (>= 1.5)")
  SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  SET(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  INCLUDE(CPack)
endif()
