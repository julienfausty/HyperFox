#ifndef HDGSOLVER_H
#define HDGSOLVER_H

#include "Solver.h"
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
  protected:
    /*!
     * \brief calculate the sparsity pattern of the matrix based on the mesh connectivity
     */
    void calcSparsityPattern();
    /*!
     * \brief a pointer to a field holding the values of the face to bulk operator for the solution
     */
    Field * U = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk operator for the flux
     */
    Field * Q = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk offset for the solution
     */
    Field * U0 = NULL;
    /*!
     * \brief a pointer to a field holding the values of the face to bulk offset for the flux
     */
    Field * Q0 = NULL;
};//HDGSolver

}//hfox

#endif//HDGSOLVER_H
