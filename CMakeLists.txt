# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(ctseg CXX Fortran)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

if (NOT ${TRIQS_WITH_PYTHON_SUPPORT})
 MESSAGE(FATAL_ERROR "Wien2TRIQS require Python support in TRIQS")
endif()

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

add_subdirectory(fortran/dmftproj) 
add_subdirectory(fortran/F90)
add_subdirectory(python)
add_subdirectory(test)
