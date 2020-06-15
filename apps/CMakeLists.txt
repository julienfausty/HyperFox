#Adds the executables.

#main executable
FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

include_directories(${MPI_CXX_INCLUDE_DIRS} ${MOAB_INCLUDE_DIRS} ${PETSC_INCLUDES} ${HDF5_C_INCLUDE_DIRS})

add_executable(
  hyperFoxDriver
  main.cpp
  )
target_link_libraries(hyperFoxDriver ${hyperFoxLibs} 
  MPI::MPI_CXX 
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  Boost::filesystem
  )


set(hyperFoxCPPFiles "")
set(hyperFoxHFiles "")
FOREACH(dir ${hyperFoxIncludes})
  file(GLOB files ${dir}/*.cpp)
  set(hyperFoxCPPFiles ${hyperFoxCPPFiles} ${files})
  file(GLOB files ${dir}/*.h)
  set(hyperFoxHFiles ${hyperFoxHFiles} ${files})
ENDFOREACH()

add_library(
  hyperfox SHARED 
  ${hyperFoxCPPFiles}
  )

#target_include_directories(hyperfox PUBLIC ${hyperFoxIncludes})

target_link_libraries(hyperfox 
  MPI::MPI_CXX 
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  Boost::filesystem
  )

install(TARGETS hyperfox EXPORT HyperFox
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/hyperFox
)

install(FILES ${hyperFoxHFiles}
  DESTINATION include/hyperFox
)

link_directories(${MOAB_LIBRARY_DIRS})

add_executable(
  tests_hyperFox
  test.cpp
  ${testSrcFiles}
  )
target_link_libraries(tests_hyperFox ${hyperFoxLibs}
  MPI::MPI_CXX
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  Boost::filesystem
  Catch2::Catch2)

add_executable(
  par_tests_hyperFox
  partest.cpp
  ${parTestSrcFiles}
  )
target_link_libraries(par_tests_hyperFox ${hyperFoxLibs}
  MPI::MPI_CXX
  ${PETSC_LIBRARIES} 
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  Boost::filesystem
  Catch2::Catch2)
