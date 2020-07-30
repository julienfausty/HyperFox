#ifndef NONLINEAROPERATOR_H
#define NONLINEAROPERATOR_H

#include "Operator.h"

namespace hfox{

/*!
 * \brief an object representing a linearized non linear operator on a field
 *
 */

class NonLinearOperator : public Operator{
  public:
    //default constructor
    using Operator::Operator;
    /*!
     * \brief get the RHS of the operator
     */
    const EVector * getRHS() const{return &rhs;};
    /*!
     * \brief a method for allocating the operator
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    virtual void allocate(int nDOFsPerNodeUser);
    /*!
     * \brief a virtual method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV the values of the determinants of jacobian matrices multiplied by the weights at the integration points
     */
    virtual void assemble(const std::vector< double > & dV, 
        const std::vector< EMatrix > & invJacobians) = 0;
  protected:
    /*!
     * \brief the rhs possible generated
     */
    EVector rhs;

};//NonLinearOperator

}//hfox

#endif//NONLINEAROPERATOR_H
