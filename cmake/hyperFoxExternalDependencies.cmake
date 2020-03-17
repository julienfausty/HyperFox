list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(MPI REQUIRED)

if (MPI_FOUND)
  set (CMAKE_C_COMPILER "${MPI_C_COMPILER}")
  set (CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
  add_definitions(-DOMPI_SKIP_MPICXX)
endif()

find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(PETSc COMPONENTS CXX)

include_directories(${PETSC_INCLUDES})

find_package(Boost REQUIRED COMPONENTS filesystem)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Eigen3 REQUIRED)

find_package(HDF5 REQUIRED)

find_package(Catch2 REQUIRED)
