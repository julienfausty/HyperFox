file(GLOB
  files
  *.h
  *.cpp
  )

add_library(parallel SHARED ${files})

target_include_directories(parallel PUBLIC ${hyperFoxIncludes} ${MPI_INCLUDE_DIRS})

target_link_libraries(parallel mesh field MPI::MPI_CXX)
