#Adds the executables.

#main executable
FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

link_directories(${MOAB_LIBRARY_DIRS})

add_executable(
  hyperFox
  main.cpp
  )
target_link_libraries(hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen ${MOAB_LIBRARIES})

add_executable(
  tests_hyperFox
  test.cpp
  ${testSrcFiles}
  )
target_link_libraries(tests_hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen ${MOAB_LIBRARIES} Catch2::Catch2)
