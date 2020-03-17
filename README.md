///////////////////////////////////////////////////////////////////////////////

                        |\                      /|
                        |  \                  /  |
                        | -  \/\/\/\/\/\/\/\/  - |
                        \        _      _        /
                          \     \-/    \./     /
                            \                /
                              \   |.T.|    /
                                \        /
                                  \    /
                                    \/

///////////////////////////////////////////////////////////////////////////////

# HyperFox Project

The HyperFox project is an attempt to create a general purpose finite element 
library for running simulations in arbitrary dimensions of both parameters and 
solutions and of arbitrary order.

# Dependencies

For now the project has the following dependencies:

- OpenMpi: For MPI parallelism support

- OpenMP: A compiler with OpenMP support should be used

- Boost: a group of libraries for tons of algorithm and container support.

- Eigen: A header library for linear algebra

- Petsc: For parallel linear solvers for linear algebra on spare matrices (with the requisite BLAS/LAPACK installations)

- HDF5: for input/output binary file formatting

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
