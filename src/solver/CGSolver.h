#ifndef CGSOLVER_H
#define CGSOLVER_H

#include "Solver.h"
#include "ErrorHandle.h"
#include "ProgressBar.h"

namespace hfox{

/*!
 * \brief The classic continuous Galerkin solver for finite elements
 */

class CGSolver : public Solver{
  public:
    /*!
     * \brief Solver constructor inheritance
     */
    using Solver::Solver;
    /*!
     * \brief a method for allocating the memory for all components
     */
    void allocate();
    /*!
     * \brief a method to assemble the linear system
     */
    void assemble();
    /*!
     * \brief a method to solve the linear system
     */
    void solve();
  protected:
    /*!
     * \brief calculate the sparsity pattern of the matrix based on the mesh connectivity
     */
    void calcSparsityPattern();
};//CGSolver

}//hfox

#endif//CGSOLVER_H
