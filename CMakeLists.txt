cmake_minimum_required(VERSION 3.4)

#Declaration of the project.

project(HyperFox)

#Project configuration.

#Variables
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Macros
include(hyperFoxMacros)

#External dependencies.
include(hyperFoxExternalDependencies)

#Add subdirectories for libraries.
add_subdirectory (src)
add_subdirectory(tests)

#Define dependencies
include(hyperFoxDependencies)

#Add subdirectory for executables.
add_subdirectory(bin)
