# Linking and include info
add_library(atm_c dos_tetra3d.hpp dos_tetra3d.cpp argsort.hpp argsort.cpp)
set_target_properties(atm_c PROPERTIES LINKER_LANGUAGE CXX)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/c++/plovasp/atm ${TRIQS_INCLUDE_ALL})
target_link_libraries(atm_c triqs)

install(TARGETS atm_c DESTINATION lib)

add_subdirectory(test)
