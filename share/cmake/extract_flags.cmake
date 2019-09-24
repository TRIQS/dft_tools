# Recursively fetch all targets that the interface of a target depends upon
macro(get_all_interface_targets name target)
  get_property(TARGET_LINK_LIBRARIES TARGET ${target} PROPERTY INTERFACE_LINK_LIBRARIES)
  foreach(lib IN LISTS TARGET_LINK_LIBRARIES)
    if(TARGET ${lib})
      # Append to list
      list(APPEND ${name}_INTERFACE_TARGETS ${lib})
      # Recure into target dependencies
      get_all_interface_targets(${name} ${lib})
    endif()
  endforeach()
endmacro()

# Extract the property from the target and recursively from all targets it depends upon
macro(get_property_recursive)
  cmake_parse_arguments(get_property_recursive "" "TARGET" "PROPERTY" ${ARGN})
  set(target ${get_property_recursive_TARGET})
  set(property ${get_property_recursive_PROPERTY})
  get_all_interface_targets(${target} ${target})
  foreach(t IN LISTS ${target}_INTERFACE_TARGETS ITEMS ${target})
    get_property(p TARGET ${t} PROPERTY ${property})
    list(APPEND ${ARGV0} ${p})
  endforeach()
  # Clean duplicates and any occurance of '/usr/include' dirs
  if(${ARGV0})
    list(REMOVE_DUPLICATES ${ARGV0})
    list(REMOVE_ITEM ${ARGV0} /usr/include)
  endif()
endmacro()

# Recursively fetch all compiler flags attached to the interface of a target
macro(extract_flags target)

  get_property_recursive(opts TARGET ${target} PROPERTY INTERFACE_COMPILE_OPTIONS)
  foreach(opt ${opts})
    set(${target}_LDFLAGS "${${target}_LDFLAGS} ${opt}")
    set(${target}_CXXFLAGS "${${target}_CXXFLAGS} ${opt}")
  endforeach()

  get_property_recursive(defs TARGET ${target} PROPERTY INTERFACE_COMPILE_DEFINITIONS)
  foreach(def ${defs})
    set(${target}_CXXFLAGS "${${target}_CXXFLAGS} -D${def}")
  endforeach()

  get_property_recursive(inc_dirs TARGET ${target} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  get_property_recursive(sys_inc_dirs TARGET ${target} PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
  list(REMOVE_ITEM inc_dirs ${sys_inc_dirs})
  foreach(dir ${inc_dirs})
    set(${target}_CXXFLAGS "${${target}_CXXFLAGS} -I${dir}")
  endforeach()
  foreach(dir ${sys_inc_dirs})
    set(${target}_CXXFLAGS "${${target}_CXXFLAGS} -isystem${dir}")
  endforeach()

  get_property_recursive(libs TARGET ${target} PROPERTY INTERFACE_LINK_LIBRARIES)
  foreach(lib ${libs})
    if(NOT TARGET ${lib})
      set(${target}_LDFLAGS "${${target}_LDFLAGS} ${lib}")
    endif()
  endforeach()

  # We have to replace generator expressions explicitly
  string(REGEX REPLACE "\\$<INSTALL_INTERFACE:([^ ]*)>" "\\1" ${target}_LDFLAGS "${${target}_LDFLAGS}")
  string(REGEX REPLACE "\\$<INSTALL_INTERFACE:([^ ]*)>" "\\1" ${target}_CXXFLAGS "${${target}_CXXFLAGS}")
  string(REGEX REPLACE " [^ ]*\\$<[^ ]*:[^ ]*>" "" ${target}_LDFLAGS "${${target}_LDFLAGS}")
  string(REGEX REPLACE " [^ ]*\\$<[^ ]*:[^ ]*>" "" ${target}_CXXFLAGS "${${target}_CXXFLAGS}")
endmacro()
