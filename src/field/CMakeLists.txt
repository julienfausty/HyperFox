file(GLOB
  files
  *.h
  *.cpp
  )

add_library(field SHARED ${files})

target_include_directories(field PUBLIC ${hyperFoxIncludes})

target_link_libraries(field mesh operator)
