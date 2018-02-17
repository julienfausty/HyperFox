file(GLOB
  files
  *.h
  *.cpp
  )

add_library(operation SHARED ${files})
