file(GLOB
  files
  *.h
  *.cpp
  )

add_library(globals SHARED ${files})

target_include_directories(globals PUBLIC ${hyperFoxIncludes} ${MPI_INCLUDE_DIRS})

target_link_libraries(globals Eigen3::Eigen MPI::MPI_CXX)
