configure_file(configRessources.h.in ${INCLUDES_OUTPUT_PATH}/configRessources.h)

SUBDIRLIST(dir_list ${CMAKE_CURRENT_SOURCE_DIR})

foreach(dir ${dir_list})
  add_subdirectory(${dir})
endforeach()
