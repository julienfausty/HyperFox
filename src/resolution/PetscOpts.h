#ifndef PETSCOPTS_H
#define PETSCOPTS_H

#include "petscksp.h"

namespace hfox{

struct PetscOpts{
  /*!
   * \brief the type of solver to be used by PetsC
   */
  KSPType solverType{ KSPGMRES };
  /*!
   * \brief the type of pre-conditionner to be used by PetsC
   */
  PCType preconditionnerType{ PCJACOBI };
  /*!
   * \brief the relative allowed tolerance for convergence to be obtained (|| rhs - mat*u || < rtol * || b ||)
   */
  PetscReal rtol{ 1e-6 };
  /*!
   * \brief the maximum number of iterations to obtain convergence
   */
  PetscInt maxits{ 1000 };
  /*!
   * \brief the verbosity of the Petsc resolution
   */
  bool verbose{ 1 };
};//PetscOpts

}//hfox

#endif //PETSCOPTS_H
