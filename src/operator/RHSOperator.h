#ifndef RHSOPERATOR_H
#define RHSOPERATOR_H

#include "Operator.h"

/*!
 * \brief an intermediate interface class specified to the right hand side
 */

namespace hfox{

class RHSOperator : public Operator{
  public:
    /*!
     * \brief constructor inheritance
     */
    using Operator::Operator;
    /*!
     * \brief allocate the operator matrix
     */
    virtual void allocate(int nDOFsPerNode);
    /*!
     * \brief a virtual method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param detJacobians the values of the determinants of jacobian matrices at the integration points
     */
    virtual void assemble(const std::vector< double > & detJacobians, 
        const std::vector< EMatrix > & invJacobians) = 0;
};//RHSOperator

}//hfox

#endif//RHSOPERATOR_H
