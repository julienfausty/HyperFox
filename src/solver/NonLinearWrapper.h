#ifndef NONLINEARWRAPPER_H
#define NONLINEARWRAPPER_H

#include <functional>
#include "Solver.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief A class that wraps around a function to iteratively solve non linear problems
 */

class NonLinearWrapper{

  public:
    /*!
     * \brief default constructor
     */
    NonLinearWrapper();
    /*!
     * \brief a method for solving the non linear problem
     */
    void solve();
    /*!
     * \brief a method for getting the residual
     */
    double getResidual(){return residual;};
    /*!
     * \brief set the residual tolerance
     * @param tol the residual tolerance being set
     */
    void setResidualTolerance(double tol){resTol = tol;};
    /*!
     * \brief set the max number of iterations
     * @param iters the maximum number of iterations
     */
    void setMaxIterations(int iters){maxIters = iters;};
    /*!
     * \brief set the current and previous solutions
     * @param currentSol the current solution at this iteration
     * @param prevSol the solution from the previous iteration
     */
    void setSolutionFields(Field * currentSol, Field * prevSol){currentSolution = currentSol; previousSolution = prevSol;};
    /*!
     * \brief set residual computer
     * @param resCom a method for computing the residual from the current and previous fields (optional)
     */
    void setResidualComputer(std::function<double(Field*, Field*)> resComp){residualComputer = resComp;};
    /*!
     * \brief set the linearized solver method
     * @param solCom a method for solving the linearized problem (mandatory)
     */
    void setLinearizedSolver(std::function<void()> solComp){linearizedSolver = solComp;};
  protected:
    /*!
     * \brief the vanilla residual computation (currentField - prevField)^2/(currentField)^2
     * @param currentSol the current solution at this iteration
     * @param prevSol the solution from the previous iteration
     */
    double vanillaResidualComputer(Field * currentSol, Field * prevSol) const;
    /*!
     * \brief a pointer to the Solver one is wrapping around
     */
    std::function<void()> linearizedSolver;
    /*!
     * \brief pointer to the previous solution
     */
    Field * previousSolution;
    /*!
     * \brief pointer to the current solution
     */
    Field * currentSolution;
    /*!
     * \brief a function to calculate the residual that takes pointers to the current solution field and previous solution field
     */
    std::function<double(Field*, Field*)> residualComputer;
    /*!
     * \brief the current residual value
     */
    double residual;
    /*!
     * \brief residual tolerance
     */
    double resTol = 1e-6;
    /*!
     * \brief a max number of iterations
     */
    int maxIters = 1000;
    /*!
     * \brief a boolean controling verbosity
     */
    bool verbosity = true;
};//NonLinearWrapper

}//hfox

#endif//NONLINEARWRAPPER_H
