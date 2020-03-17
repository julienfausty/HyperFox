#ifndef SOLVER_H
#define SOLVER_H

#include "LinAlgebraInterface.h"
#include "Model.h"

namespace hfox{

/*!
 * \brief The interface to create solvers
 *
 * The main goal of this interface to to be able to assemble and solve systems 
 * independently of the models that are being used.
 */

class Solver{
  public:
    /*!
     * \brief a method to set the model
     */
    virtual void setModel(Model * m){model = m;};
    /*!
     * \brief a method to set the linear algebra interface
     */
    virtual void setLinSystem(LinAlgebraInterface * lai){LinSystem = lai;};
    /*!
     * \brief a method to assemble the linear system
     */
    virtual void assemble()=0;
    /*!
     * \brief a method to solve the linear system
     */
    virtual void solve()=0;
  protected:
    /*!
     * \brief the interface to the linear algebra package
     */
    LinAlgebraInterface * LinSystem;
    /*!
     * \brief the model of the system being simulated
     */
    Model * model;

};//Solver

}//hfox

#endif //SOLVER_H
