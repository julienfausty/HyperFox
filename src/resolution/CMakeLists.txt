file(GLOB
  files
  *.h
  *.cpp
  )

add_library(resolution SHARED ${files})

target_include_directories(resolution PUBLIC ${hyperFoxIncludes} ${PETSC_INCLUDES} ${MPI_INCLUDE_DIRS})

target_link_libraries(resolution globals ${PETSC_LIBRARIES} MPI::MPI_CXX)
