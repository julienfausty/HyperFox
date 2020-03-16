#Adds the executables.

#main executable
FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

add_executable(
  hyperFox
  main.cpp
  )
target_link_libraries(hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen)

add_executable(
  tests_hyperFox
  test.cpp
  ${testSrcFiles}
  )
target_link_libraries(tests_hyperFox ${hyperFoxLibs} MPI::MPI_CXX ${PETSC_LIBRARIES} Eigen3::Eigen Catch2::Catch2)
