SUBDIRLIST(dir_list ${CMAKE_CURRENT_SOURCE_DIR})

FOREACH(dir ${dir_list})
  add_subdirectory(${dir})
ENDFOREACH()
