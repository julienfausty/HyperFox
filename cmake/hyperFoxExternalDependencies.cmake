add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE $ENV{CATCH_DIR}/single_include)
