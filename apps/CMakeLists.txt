#Adds the executables.

#main executable
FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

include_directories(${MPI_CXX_INCLUDE_DIRS} ${MOAB_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

add_executable(
  hyperFox
  main.cpp
  )
target_link_libraries(hyperFox ${hyperFoxLibs} 
  MPI::MPI_C 
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  )

add_executable(
  tests_hyperFox
  test.cpp
  ${testSrcFiles}
  )
target_link_libraries(tests_hyperFox ${hyperFoxLibs}
  MPI::MPI_C 
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  Catch2::Catch2)
