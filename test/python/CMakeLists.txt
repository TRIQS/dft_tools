# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(all_tests Toto chain)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# If ASAN=ON use minimal runtime of UBSAN if it exists
if(${ASAN} AND UBSAN_MINIMAL_RT_LIBRARY)
  set(UBSAN_RT_LIBRARY ${UBSAN_MINIMAL_RT_LIBRARY})
endif()

# Set the PythonPath, Sanitizer Library Preloads and UBSan Options for all tests
set_property(TEST ${all_tests} PROPERTY
  ENVIRONMENT
    PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
    LD_PRELOAD=$<$<BOOL:${ASAN}>:${ASAN_RT_LIBRARY}>:$<$<BOOL:${UBSAN}>:${UBSAN_RT_LIBRARY}>
    UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
)
