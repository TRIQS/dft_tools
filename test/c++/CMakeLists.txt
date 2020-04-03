# Build googletest
set(PYTHON_EXECUTABLE ${TRIQS_PYTHON_EXECUTABLE})
add_subdirectory(gtest EXCLUDE_FROM_ALL)

# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  target_link_libraries(${test_name} app4triqs_c gtest_main project_warnings)
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++17
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
  endif()
endforeach()
