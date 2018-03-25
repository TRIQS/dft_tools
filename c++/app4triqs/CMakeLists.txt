file(GLOB_RECURSE sources *.cpp *.hpp)
add_library(app4triqs_c ${sources})

target_link_libraries(app4triqs_c PUBLIC triqs)
target_compile_options(app4triqs_c PUBLIC -std=c++17)
target_include_directories(app4triqs_c PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/c++>)
target_compile_definitions(app4triqs_c PRIVATE
				APP4TRIQS_GIT_HASH=${APP4TRIQS_GIT_HASH}
			    	TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:APP4TRIQS_DEBUG>
			        $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			 	)

# Install library and headers
install(TARGETS app4triqs_c EXPORT app4triqs-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# ========= Static Analyzer Checks ==========

# Locate static analyzer tools
find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

# Run clang-tidy if found
if(CLANG_TIDY_EXECUTABLE)
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
  set_target_properties(app4triqs_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
else()
  message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
endif()

# Run cppcheck if found
if(CPPCHECK_EXECUTABLE)
  message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
  add_custom_command(
	TARGET app4triqs_c
	COMMAND echo "--- Running cppcheck ---\n"
	COMMAND ${CPPCHECK_EXECUTABLE}
	--enable=warning,style,performance,portability
        --std=c++14
	--template=gcc
        --verbose
        --quiet
        ${sources}
	COMMAND echo "------------------------\n"
 )
else()
  message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
endif()


# ========= Dynamic Analyzer Checks ==========

# Address Sanitizer
option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(FATAL_ERROR "Address Sanitizer Checks currently not compatible with OSX.")
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/share/cmake)
  find_package(libasan_rt REQUIRED)
  target_compile_options(app4triqs_c PUBLIC -fsanitize=address -fno-omit-frame-pointer)
  target_link_libraries(app4triqs_c INTERFACE "-fsanitize=address -fno-omit-frame-pointer")
endif()

# Undefined Behavior Sanitizer
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(FATAL_ERROR "Undefined Behavior Sanitizer Checks currently not compatible with OSX.")
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/share/cmake)
  find_package(libubsan_rt REQUIRED)
  target_compile_options(app4triqs_c PUBLIC -fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize=vptr)
  target_link_libraries(app4triqs_c INTERFACE "-fsanitize=undefined -fno-omit-frame-pointer")
endif()
