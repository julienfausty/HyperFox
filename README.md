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

For now the project has two dependencies:

- OpenMpi: For MPI parallelism support.

- OpenMP: A compiler with OpenMP support should be used.

- Eigen: A header library for linear algebra.

- Catch2: A unittesting library

In the future, the project will most likely depend on:

- Petsc: For parallel linear solvers for linear algebra on spare matrices (with the requisite BLAS/LAPACK installations)
