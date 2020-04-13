# HyperFox Project

The HyperFox project is an attempt to create a general purpose finite element 
library for running simulations in arbitrary dimensions of both parameters and 
solutions and of arbitrary order with favoring of DG methods.

# Dependencies

For now the project has the following dependencies:

- CMake: for generating the makefiles easily

- OpenMpi: For MPI parallelism support

- OpenMP: A compiler with OpenMP support should be used

- Boost: a group of libraries for tons of algorithm and container support.

- Eigen: A header library for linear algebra

- Petsc: For parallel linear solvers for linear algebra on spare matrices (with the requisite BLAS/LAPACK installations, with clanguage=cxx)

- HDF5: for input/output binary file formatting

- MOAB: a c++ library for dealing with meshes and computing mesh entities.

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

The tests are all linked to the "tests_hyperfox" binary. Tu run all the tests one can simply run

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
