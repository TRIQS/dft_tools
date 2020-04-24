include(external_dependency.cmake)

# Add your dependencies with the function
#
#   external_dependency(name
#       [VERSION <version-number>]
#       [GIT_REPO <url>]
#       [TAG <tag>]
#       [EXCLUDE_FROM_ALL]
#   )
#
# Resolve the dependency using the following steps in order.
# If a step was successful, skip the remaining ones.
#
#  1. Use find_package(name [<version-number>])
#     to locate the package in the system.
#     Skip this step if Build_Deps option is set.
#  2. Try to find a directory containing the sources
#     at ${CMAKE_SOURCE_DIR}/deps/name. If found
#     build it as a cmake sub-project.
#  3. If GIT_REPO is provided, git clone the sources,
#     and build them as a cmake sub-project.
#
# Addtional options:
#
#  TAG - Use this keyword to specify the git tag or branch
#
#  EXCLUDE_FROM_ALL - If set, targets of the dependency cmake subproject
#                     will not be included in the ALL target of the project.
#                     In particular the dependency will not be installed.

option(Build_Deps "Build dependencies from source" OFF)

# -- GTest --
external_dependency(GTest
  GIT_REPO github:google/googletest
  TAG release-1.10.0
  EXCLUDE_FROM_ALL
)
