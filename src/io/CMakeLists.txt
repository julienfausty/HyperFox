file(GLOB
  files
  *.h
  *.cpp
  )

add_library(io SHARED ${files})
