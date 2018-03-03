add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE $ENV{CATCH_DIR}/single_include)

add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE 
  $ENV{EIGEN_DIR}/include/eigen3)
