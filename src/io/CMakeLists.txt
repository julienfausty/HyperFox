file(GLOB
  files
  *.h
  *.cpp
  )

add_library(io SHARED ${files})

target_include_directories(io PUBLIC ${hyperFoxIncludes} ${MOAB_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

target_link_libraries(io globals mesh field Boost::filesystem ${MOAB_LIBRARIES} ${HDF5_C_LIBRARIES} MPI::MPI_CXX)
