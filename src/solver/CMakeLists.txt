file(GLOB
  files
  *.h
  *.cpp
  )

add_library(solver SHARED ${files})

target_include_directories(solver PUBLIC ${hyperFoxIncludes})

target_link_libraries(solver model resolution mesh field parallel)
