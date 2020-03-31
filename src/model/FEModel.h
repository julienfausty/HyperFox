#ifndef FEMODEL_H
#define FEMODEL_H

#include <map>
#include <string>
#include "Model.h"
#include "ReferenceElement.h"
#include "DenseEigen.h"
#include "Operator.h"


namespace hfox{

/*!
 * \brief An interface class that describes the necessary elements of a finite element model
 *
 * The finite element model should implement the relevant physics towards the assembly of the linear system
 * for solving the finite element problem.
 */

class FEModel : public Model{
  public:
    /*!
     * \brief the implementation of the compute method
     */
    virtual void compute();
    /*!
     * \brief compute the local matrix
     */
    virtual const EMatrix getLocalMatrix() const{return localMatrix;};
    /*!
     * \brief compute the local right hand side
     */
    virtual const EVector getLocalRHS() const{return localRHS;};
  protected:
    /*!
     * \brief a method where the relevant operators should be initialized
     */
    virtual void initializeOperators()=0;
    /*!
     * \brief compute the local matrix
     */
    virtual void computeLocalMatrix()=0;
    /*!
     * \brief compute the local right hand side
     */
    virtual void computeLocalRHS()=0;
    /*!
     * \brief the local element matrix
     */
    EMatrix localMatrix;
    /*!
     * \brief the local element right hand side
     */
    EVector localRHS;
    /*!
     * \brief a map containing the local operators relevant to the model with their names
     */
    std::map<std::string, Operator> operatorMap;

}; //FEModel

}//hfox

#endif //FEMODEL_H
