file(GLOB
  files
  *.h
  *.cpp
  )

add_library(mesh SHARED ${files})
