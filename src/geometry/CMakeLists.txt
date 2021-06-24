file(GLOB
  files
  *.h
  *.cpp
  )

add_library(geometry SHARED ${files})

target_include_directories(geometry PUBLIC ${hyperFoxIncludes})

target_link_libraries(geometry globals element operator mesh field)
