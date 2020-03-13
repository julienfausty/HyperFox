file(GLOB
  files
  *.h
  *.cpp
  )

add_library(mesh SHARED ${files})

target_include_directories(mesh PUBLIC ${hyperFoxIncludes})

target_link_libraries(mesh io Eigen3::Eigen)
