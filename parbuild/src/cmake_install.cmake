# Install script for directory: /home/jfausty/workspace/Codes/HyperFox/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/element/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/field/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/globals/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/io/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/mesh/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/model/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/operator/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/parallel/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/resolution/cmake_install.cmake")
  include("/home/jfausty/workspace/Codes/HyperFox/parbuild/src/solver/cmake_install.cmake")

endif()

