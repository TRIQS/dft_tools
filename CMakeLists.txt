# Version number
set(APP4TRIQS_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 3.6)
project(app4triqs CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Build static libraries
set(BUILD_SHARED_LIBS OFF)

# Enable compiler warnings for the whole project
add_definitions(-Wall -fPIC)

# Load Dependencies
find_package(TRIQS ${APP4TRIQS_VERSION} EXACT REQUIRED)
find_package(Cpp2Py ${APP4TRIQS_VERSION} EXACT REQUIRED)

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Fetch the git hash of app4triqs
triqs_get_git_hash_of_source_dir(APP4TRIQS_GIT_HASH)
message(STATUS "app4triqs git hash: ${APP4TRIQS_GIT_HASH}")

# Build and install the app4triqs library
add_subdirectory(c++/app4triqs)

# Build and install the app4triqs python module
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
