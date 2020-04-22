# Each dependency should be either
#  * a cmake subproject in this directory
#  * found in the system using find_package
# Provide either a download.sh script to get all sources
# or use e.g. git subtree/submodule or softlinks

# Add all subdirectories as CMake Sub-Projects
file(GLOB allfiles CONFIGURE_DEPENDS *)
foreach(file ${allfiles})
  if(IS_DIRECTORY ${file})
    add_subdirectory(${file})
  endif()
endforeach()

# -- googletest --
if(NOT TARGET gtest)
  find_package(GTest REQUIRED)
endif()
