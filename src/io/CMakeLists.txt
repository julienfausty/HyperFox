file(GLOB
  files
  *.h
  *.cpp
  )

add_library(io SHARED ${files})

target_include_directories(io PUBLIC ${hyperFoxIncludes} ${MOAB_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

target_link_libraries(io mesh field Boost::filesystem ${MOAB_LIBRARIES})
