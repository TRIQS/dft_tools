# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(all_tests Toto chain)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath and the Sanitizer Library Preloads
set_property(TEST ${all_tests} PROPERTY
  ENVIRONMENT
    PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
    LD_PRELOAD=$<$<BOOL:${ASAN}>:${ASAN_RT_LIBRARY}>:$<$<BOOL:${UBSAN}>:${UBSAN_RT_LIBRARY}>
)
