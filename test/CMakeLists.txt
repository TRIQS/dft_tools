add_subdirectory(c++)

if(PythonSupport)
 add_subdirectory(python)
 add_subdirectory(python/plovasp)
 add_subdirectory(python/elk)
endif()
