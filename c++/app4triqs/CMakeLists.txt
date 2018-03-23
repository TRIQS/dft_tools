file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(app4triqs_c ${sources})

target_link_libraries(app4triqs_c PUBLIC triqs)
target_compile_options(app4triqs_c PUBLIC -std=c++17)
target_include_directories(app4triqs_c PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/c++>)

target_compile_definitions(app4triqs_c PRIVATE
				APP4TRIQS_GIT_HASH=${APP4TRIQS_GIT_HASH}
			    	TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:DEBUG_APP4TRIQS>
			        $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
				#$<$<BOOL:${MyOption}>:MYOPTION> # Add compile definintions for option
			 	)

# Install
install(TARGETS app4triqs_c EXPORT app4triqs-targets DESTINATION lib)

# If other applications need our headers, we can install them as well
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
