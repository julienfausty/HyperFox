file(GLOB
  files
  *.h
  *.cpp
  )

add_library(mesh SHARED ${files})

target_include_directories(mesh PUBLIC ${hyperFoxIncludes} ${MOAB_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

target_link_libraries(mesh element Eigen3::Eigen ${MOAB_LIBRARIES} MPI::MPI_CXX)
