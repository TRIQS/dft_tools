include(get_cpm.cmake)

# -- dftkit --
CPMAddPackage(
  NAME dftkit
  GITHUB_REPOSITORY triqs/dftkit
  GIT_TAG unstable
)
