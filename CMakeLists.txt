cmake_minimum_required(VERSION 3.4)

#Declaration of the project.

project(HyperFox LANGUAGES CXX C VERSION 0.0.1)

#Project configuration.

set (CMAKE_CXX_STANDARD 17)

set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting default build type: Release")
  set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE} CACHE STRING
          "Choose the type of build, options are: Debug Release
          RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)

#cmake config
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")

#Variables
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(INCLUDES_OUTPUT_PATH ${PROJECT_BINARY_DIR}/includes)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Macros
include(hyperFoxMacros)

#External dependencies.
include(hyperFoxExternalDependencies)

#Define dependencies
include(hyperFoxDependencies)

#Add subdirectories for libraries.
add_subdirectory (src)
add_subdirectory(tests)

#Add subdirectory for executables.
add_subdirectory(apps)

#Add subdirectory for tools.
add_subdirectory(tools)

# Add all targets to the build-tree export set
export(TARGETS hyperfox
  FILE "${PROJECT_BINARY_DIR}/HyperFox.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE HyperFox)

configure_file(cmake/HyperFoxConfig.cmake.in
  "${PROJECT_BINARY_DIR}/HyperFoxConfig.cmake" @ONLY)

configure_file(cmake/HyperFoxConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/HyperFoxConfigVersion.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/HyperFoxConfig.cmake"
  "${PROJECT_BINARY_DIR}/HyperFoxConfigVersion.cmake"
  DESTINATION share)

# Install the export set for use with the install-tree
install(EXPORT HyperFox DESTINATION
  "${CMAKE_INSTALL_PREFIX}/share")


