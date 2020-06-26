# Install script for directory: /home/jfausty/workspace/Codes/HyperFox/apps

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/jfausty/workspace/Codes/HyperFox/parbuild/lib/libhyperfox.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so"
         OLD_RPATH "/home/jfausty/workspace/Dependencies/petsc-3.13.2/arch-linux-c-debug/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/local/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhyperfox.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/hyperFox" TYPE FILE FILES
    "/home/jfausty/workspace/Codes/HyperFox/src/element/Cubature.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/element/ElementGeometry.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/element/ReferenceElement.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/field/Field.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/field/FieldTypes.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/globals/DenseEigen.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/globals/ErrorHandle.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/globals/Modifier.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/globals/ProgressBar.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/globals/Utils.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/io/HDF5Io.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/io/Io.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/io/MoabMeshIo.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/mesh/Mesh.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/AssemblyType.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/DiffusionSource.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/DirichletModel.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/FEModel.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/HDGConvectionDiffusionReactionSource.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/HDGDiffusionSource.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/HDGLaplaceModel.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/HDGModel.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/HDGTransport.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/LaplaceModel.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/Model.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/model/Transport.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Convection.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Diffusion.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Euler.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/HDGBase.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/HDGConvection.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/HDGDiffusion.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/HDGOperator.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Mass.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Operator.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/RHSOperator.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/RKType.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Reaction.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/RungeKutta.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/Source.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/operator/TimeScheme.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/parallel/Partitioner.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/parallel/ZoltanOpts.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/parallel/ZoltanPartitioner.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/resolution/LinAlgebraInterface.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/resolution/PetscInterface.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/resolution/PetscOpts.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/solver/CGSolver.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/solver/HDGSolver.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/solver/HDGSolverOpts.h"
    "/home/jfausty/workspace/Codes/HyperFox/src/solver/Solver.h"
    "/home/jfausty/workspace/Codes/HyperFox/parbuild/includes/TestUtils.h"
    )
endif()

