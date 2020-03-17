file(GLOB
  files
  *.h
  *.cpp
  )

add_library(operation SHARED ${files})

target_include_directories(operation PUBLIC ${hyperFoxIncludes})

target_link_libraries(operation field)
