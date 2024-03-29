#ifndef HDGSOLVER_H
#define HDGSOLVER_H

#include <algorithm>
#include "Solver.h"
#include "HDGSolverOpts.h"
#include "ErrorHandle.h"
#include "ProgressBar.h"

namespace hfox{

/*!
 * \brief General solver for Hybrid Discontinuous Galerkin finite elements.
 */

class HDGSolver : public Solver{
  public:
    /*!
     * \brief Solver constructor inheritance
     */
    using Solver::Solver;
    /*!
     * \brief destructor
     */
    ~HDGSolver();
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
    /*!
     * \brief a method to set the options of the solver
     */
    void setOptions(HDGSolverOpts opts);
  protected:
    /*!
     * \brief calculate the sparsity pattern of the matrix based on the mesh connectivity
     */
    void calcSparsityPattern();
    /*!
     * \brief a method to calculate the elemental matrices
     */
    void calcElementalMatrices();
    /*!
     * \brief a method that applies boundary conditions to the elemental matrices
     */
    void applyBoundaryConditions();
    /*!
     * \brief assemble the linear system from S and S0
     */
    void assembleSystem();
    /*!
     * \brief a pointer to a field holding the values of the face to bulk operator for the solution
     */
    Field * U = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk operator for the flux
     */
    Field * Q = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to face operator for the trace
     */
    Field * S = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk offset for the solution
     */
    Field * U0 = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk offset for the flux
     */
    Field * Q0 = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to face offset for the trace
     */
    Field * S0 = NULL;
    /*!
     * \brief structure for setting options of solver
     */
    HDGSolverOpts myOpts;
};//HDGSolver

}//hfox

#endif//HDGSOLVER_H
