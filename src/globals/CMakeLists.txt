file(GLOB
  files
  *.h
  *.cpp
  )

add_library(globals SHARED ${files})

target_include_directories(globals PUBLIC ${hyperFoxIncludes})

target_link_libraries(globals)
