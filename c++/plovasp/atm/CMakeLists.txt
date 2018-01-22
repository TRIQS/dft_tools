add_library(atm_c dos_tetra3d.hpp dos_tetra3d.cpp argsort.hpp argsort.cpp)
target_link_libraries(atm_c triqs)
target_compile_options(atm_c PRIVATE -std=c++17)

install(TARGETS atm_c DESTINATION lib)

add_subdirectory(test)
