# Copy reference h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(all_tests toto)

foreach(t ${all_tests})
  add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  target_link_libraries(${t} app4triqs_c gtest)
  add_test(${t} ${CMAKE_CURRENT_BINARY_DIR}/${t})
endforeach()
