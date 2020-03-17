file(GLOB
  files
  *.h
  *.cpp
  )

add_library(io SHARED ${files})

target_include_directories(io PUBLIC ${hyperFoxIncludes})

target_link_libraries(io mesh field Boost::filesystem)
