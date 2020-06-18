# HyperFox Project

The HyperFox project is an attempt to create a general purpose finite element 
library for running simulations in arbitrary dimensions of both parameters and 
solutions and of arbitrary order favoring of Discontinuous Garlerkin (DG) methods.

The project is divided into multiple interdependent libraries.

The mesh is encapsulated by the Mesh class. Fields, defined on the mesh, are encapsulated by 
the field class which is flexible enough to describe continuous and discontinuous fields.

The Solver class implements Solver objects that assemble the global linear system from a 
Model object. The Model objects are meant to represent the finite element system of 
equations in an abstract way and define the local element by element assembly. The Model 
class uses Operators as building blocks to implement the most common finite element terms
(i.e Diffusion, Convection, Source, etc.). The general reference element of arbitrary 
dimension, order and geometry is implemented in the ReferenceElement class.

Once a system is assembled, it is solved through the use of a LinAlgebraInterface (owned by 
the Solver) which is an interface class to an external linear algebra package.

Input of meshes and parameters as well as output of solution fields are taken care of by the 
Io classes.

The idea behind the structure of the software is to sufficiently separate all the unitary 
operations and achieve a certain level of abstraction. As such, the development of a new 
Model does not need the rewritting of a Solver. Including a new element geometry 
should not generate any issues in the already implemented classes. The parallelism of the 
code only comes into play at the Solver and Mesh levels and thus writing new models for existing 
Solvers can be done rather transparently in serial. The entire software should remain modular.

In a context where developpers are usually active for only three years and at a high turnover rate, 
the structure of the code should ensure its longevity.

Every new implementation should have its methods unittested in the "tests/unittests" directory. 
Also, when applicable, small test cases should be included in the "tests/regression" directory 
for ensuring higher level behavior of the code. This test database garantees the robustness of the 
code by failing whenever an aspect of the structure is compromised. The "git" versionning software 
allows for mistakes to be made and unmade at no cost to the integrity of the application.

# Goals

- Short term (~ weeks):
  - Parallelize the code: the plan is to link to an external partitionner (Zoltan or Metis) to partition the 
  mesh and then to adapt the Solvers to parallel assembly (the current linear algebra package PetsC is 
  already parallel)
  - Re-visit the previous explicit/implicit HDG study in parallel

- Long term (~ months):
  - Implement plasma models
  - Adapt parallelism for more heterogeneous architectures (maybe GPU acceleration? PetsC already integrates 
  cuda and OpenCL specific code -> towards intensive computing for 3D)

# Dependencies

For now the project has the following dependencies:

- CMake: for generating the makefiles easily

- OpenMpi: For MPI parallelism support (check that mpi-io capabilities are present for full parallel hdf5 support)

- OpenMP: A compiler with OpenMP support should be used

- Boost: a group of libraries for tons of algorithm and container support.

- Eigen: A header library for linear algebra

- Petsc: For iterative parallel linear solvers for linear algebra on spare matrices (with the requisite BLAS/LAPACK installations, with clanguage=cxx)

- HDF5: for input/output binary file formatting

- MOAB: a c++ library for dealing with meshes and computing mesh entities 

- Zoltan: a c++ library for graph partitioning (must set ZOLTAN_PREFIX envioronment variable to install directory of Zoltan)

Also, the following library is used for the testing module:

- Catch2: A unittesting library

All these dependencies are listed as required and use the 'find_package' utility given by cmake.

# Making the project

The project is built with cmake so the following commands should create the libraries and binaries, if all the requirements are met:

```bash
mkdir build
cd build
cmake ..
make
```

# To run tests

The tests are all linked to the "tests_hyperfox" binary. To run all the tests one can simply run

```bash
./bin/tests_hyperfox
```

from the build directory. To run tests specific to one library of the project run

```bash
./bin/tests_hyperfox [library_name]
```

In order to run tests for one specific class run

```bash
./bin/tests_hyperfox [class_name]
```
