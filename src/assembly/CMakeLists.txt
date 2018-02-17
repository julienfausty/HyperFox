file(GLOB
  files
  *.h
  *.cpp
  )

add_library(assembly SHARED ${files})
