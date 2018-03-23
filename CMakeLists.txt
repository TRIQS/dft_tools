# Version number
set(APP4TRIQS_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 2.8)
project(app4triqs CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS ${APP4TRIQS_VERSION} EXACT REQUIRED)
find_package(Cpp2Py ${APP4TRIQS_VERSION} EXACT REQUIRED)

# We will need the FindNNFT of triqs, so we adjust CMake Module Path to find it
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# Default Install directory to TRIQS_ROOT if not given. Checks an absolute name is given.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Fetch the git hash of app4triqs
triqs_get_git_hash_of_source_dir(APP4TRIQS_GIT_HASH)
message(STATUS "app4triqs git hash: ${APP4TRIQS_GIT_HASH}")

# Set up the rpath for compiled executable
triqs_set_rpath_for_all_targets()

add_subdirectory(c++/app4triqs)
if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python/app4triqs)
endif()

# Tests
enable_testing()
add_subdirectory(test)

# Additional configuration files
add_subdirectory(share)

option(Build_Documentation "Build documentation" OFF)
if(${Build_Documentation})
  if(NOT ${TRIQS_WITH_DOCUMENTATION})
    message("Error: TRIQS library has not been compiled with its documentation")
  else()
    add_subdirectory(doc)
  endif()
endif()
