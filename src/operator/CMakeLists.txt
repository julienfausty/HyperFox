file(GLOB
  files
  *.h
  *.cpp
  )

add_library(operator SHARED ${files})

target_include_directories(operator PUBLIC ${hyperFoxIncludes})

target_link_libraries(operator globals element)
