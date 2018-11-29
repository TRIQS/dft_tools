# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(all_tests Toto chain)

# Nasty Workaround for Mac Os Homebrew Python to find the proper Sanitizer Compatible Python Executable
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND (ASAN OR UBSAN))
  execute_process(COMMAND bash "-c" "readlink -f $(which python) | xargs dirname" OUTPUT_VARIABLE PYTHON_EXECUTABLE OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PYTHON_EXECUTABLE "${PYTHON_EXECUTABLE}/../Resources/Python.app/Contents/MacOS/Python")
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    message(FATAL_ERROR "Could not find the proper Python executable for MacOs")
  endif()
else()
  set(PYTHON_EXECUTABLE python)
endif()

message(STATUS "Python Executable: ${PYTHON_EXECUTABLE}")

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})

# === Preload Sanitizer Libraries ===

# If ASAN=ON use minimal runtime of UBSAN if it exists
if(ASAN AND UBSAN_MINIMAL_RT_LIBRARY)
  set(UBSAN_RT_LIBRARY ${UBSAN_MINIMAL_RT_LIBRARY})
endif()

# Set the Sanitizer Library Preloads and UBSan Options for all tests
if(ASAN OR UBSAN)
  set_property(TEST ${all_tests} PROPERTY
    ENVIRONMENT
      PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
      LD_PRELOAD=$<$<BOOL:${ASAN}>:${ASAN_RT_LIBRARY}>:$<$<BOOL:${UBSAN}>:${UBSAN_RT_LIBRARY}>
      DYLD_INSERT_LIBRARIES=$<$<BOOL:${ASAN}>:${ASAN_RT_LIBRARY}>:$<$<BOOL:${UBSAN}>:${UBSAN_RT_LIBRARY}>
      UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1:$ENV{UBSAN_OPTIONS}
  )
endif()
