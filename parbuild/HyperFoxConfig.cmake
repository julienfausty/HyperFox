# - Config file for the FooBar package
# It defines the following variables
#  HYPERFOX_INCLUDES - include directories for HyperFox
#  HYPERFOX_LIBRARIES    - libraries to link against

include(CMakeFindDependencyMacro)

find_dependency(MPI REQUIRED)

find_dependency(OpenMP REQUIRED)

find_dependency(PETSc COMPONENTS CXX)

find_dependency(Boost REQUIRED COMPONENTS filesystem)

find_dependency(Eigen3 REQUIRED)

find_dependency(HDF5 REQUIRED)

find_dependency(Catch2 REQUIRED)

find_dependency(MOAB REQUIRED)

find_dependency(Zoltan REQUIRED)

# Compute paths
get_filename_component(HYPERFOX_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(HYPERFOX_INCLUDES "/usr/local/include")

# Our library dependencies (contains definitions for IMPORTED targets)

include("${CMAKE_CURRENT_LIST_DIR}/HyperFox.cmake")


# These are IMPORTED targets created by HyperFoxTargets.cmake
set(HYPERFOX_LIBRARIES hyperfox)

