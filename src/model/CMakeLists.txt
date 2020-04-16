file(GLOB
  files
  *.h
  *.cpp
  )

add_library(model SHARED ${files})

target_include_directories(model PUBLIC ${hyperFoxIncludes})

target_link_libraries(model globals element operator)
