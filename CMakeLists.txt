cmake_minimum_required(VERSION 3.4)

#Declaration of the project.

project(HyperFox LANGUAGES CXX C VERSION 0.0.1)

#Project configuration.

set (CMAKE_CXX_STANDARD 17)

#Variables
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(INCLUDES_OUTPUT_PATH ${PROJECT_BINARY_DIR}/includes)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Macros
include(hyperFoxMacros)

#External dependencies.
include(hyperFoxExternalDependencies)

set (CMAKE_C_COMPILER "${MPI_C_COMPILER}")
set (CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

#Define dependencies
include(hyperFoxDependencies)

#Add subdirectories for libraries.
add_subdirectory (src)
add_subdirectory(tests)

#Add subdirectory for executables.
add_subdirectory(apps)
