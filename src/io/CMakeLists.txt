file(GLOB
  files
  *.h
  *.cpp
  )

add_library(io SHARED ${files})

target_link_libraries(io Boost::filesystem)
