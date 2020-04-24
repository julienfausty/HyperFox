#ifndef HDGOPERATOR_H
#define HDGOPERATOR_H

#include "Operator.h"

namespace hfox{

/*!
 * \brief Interface class to Hybrid Discontinuous Galerkin type operators
 */

class HDGOperator : public Operator{

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
     * @param dV values of the discrete measure
     */
    virtual void assemble(const std::vector< double > & dV, 
        const std::vector< EMatrix > & invJacobians) = 0;

};//HDGOperator

}//hfox


#endif //HDGOPERATOR_H
