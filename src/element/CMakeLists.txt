file(GLOB
  files
  *.h
  *.cpp
  )

add_library(element SHARED ${files})

target_include_directories(element PUBLIC ${hyperFoxIncludes})

target_link_libraries(element globals)
