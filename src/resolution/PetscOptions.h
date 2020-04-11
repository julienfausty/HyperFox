#ifndef PETSCOPTIONS_H
#define PETSCOPTIONS_H

#include "petscksp.h"

namespace hfox{

struct PetscOptions{
  /*!
   * \brief the type of solver to be used by PetsC
   */
  KSPType solverType = KSPGMRES;
  /*!
   * \brief the type of pre-conditionner to be used by PetsC
   */
  PCType preconditionnerType = PCJACOBI;
  /*!
   * \brief the relative allowed tolerance for convergence to be obtained (|| rhs - mat*u || < rtol * || b ||)
   */
  PetscReal rtol = 1e-6;
  /*!
   * \brief the maximum number of iterations to obtain convergence
   */
  PetscInt maxits = 1e3;
  /*!
   * \brief the verbosity of the Petsc resolution
   */
  bool verbose = 1;
};//PetscOptions

}//hfox

#endif //PETSCOPTIONS_H
