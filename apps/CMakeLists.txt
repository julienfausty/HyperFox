#Adds the executables.

#main executable
FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

include_directories(${MPI_CXX_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

add_executable(
  hyperFox
  main.cpp
  )
target_link_libraries(hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen ${MOAB_LIBRARIES} MPI::MPI_CXX)

add_executable(
  tests_hyperFox
  test.cpp
  ${testSrcFiles}
  )
target_link_libraries(tests_hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen ${MOAB_LIBRARIES} MPI::MPI_CXX Catch2::Catch2)
