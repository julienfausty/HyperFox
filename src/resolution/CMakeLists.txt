file(GLOB
  files
  *.h
  *.cpp
  )

add_library(resolution SHARED ${files})
