#ifndef NONLINEARWRAPPER_H
#define NONLINEARWRAPPER_H

#include <functional>
#include <mpi.h>
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
     * \brief destructor
     */
    ~NonLinearWrapper();
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
     * \brief set the verbosity
     * @param verbosity boolean for verbose or not
     */
    void setVerbosity(bool verbosity){verbose = verbosity;};
    /*!
     * \brief set the current and previous solutions
     * @param currentSol the current solution at this iteration
     * @param prevSol the solution from the previous iteration
     */
    void setSolutionFields(Field * currentSol, Field * prevSol){currentSolution = currentSol; previousSolution = prevSol;};
    /*!
     * \brief set the solver
     * @param solver the solver
     */
    void setSolver(Solver * solver){mySolver = solver;};
    /*!
     * \brief set residual computer
     * @param resCom a method for computing the residual from the current and previous fields (optional)
     */
    void setResidualComputer(std::function<double(Field*, Field*)> resComp){residualComputer = resComp;};
    /*!
     * \brief set the linearized solver method
     * @param solCom a method for solving the linearized problem (optional)
     */
    void setLinearizedSolver(std::function<void(Solver*)> solComp){linearizedSolver = solComp;};
    /*!
     * \brief method to set dampening coefficient
     * @param damp the suggested value
     */
    void setDampening(double damp){dampening = damp;};
  protected:
    /*!
     * \brief the vanilla residual computation (currentField - prevField)^2/(currentField)^2
     * @param currentSol the current solution at this iteration
     * @param prevSol the solution from the previous iteration
     */
    static double vanillaResidualComputer(Field * currentSol, Field * prevSol);
    /*!
     * \brief the vanilla solve
     * @param solver the solver
     */
    static void vanillaLinearizedSolver(Solver * solver);
    /*!
     * \brief a pointer to the Solver one is wrapping around
     */
    std::function<void(Solver*)> linearizedSolver;
    /*!
     * \brief the solver to solve for
     */
    Solver * mySolver;
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
    int maxIters = 20;
    /*!
     * \brief a boolean controling verbosity
     */
    bool verbose = true;
    /*!
     * \brief a positive value used to dampen or overconverge the solution
     */
    double dampening = 0.0;
};//NonLinearWrapper

}//hfox

#endif//NONLINEARWRAPPER_H
