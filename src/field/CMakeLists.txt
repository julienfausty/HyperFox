file(GLOB
  files
  *.h
  *.cpp
  )

add_library(field SHARED ${files})
