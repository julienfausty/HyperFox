FOREACH(dir ${hyperFoxIncludes})
  include_directories(${dir})
ENDFOREACH()

include_directories(${MOAB_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS})

link_directories(${MOAB_LIBRARY_DIRS})

add_executable(
  convertGmsh2H5HO
  convertGmsh2H5HO.cpp
  )

target_link_libraries(convertGmsh2H5HO ${hyperFoxLibs}
  Eigen3::Eigen 
  ${MOAB_LIBRARIES} 
  ${HDF5_C_LIBRARIES} 
  )
