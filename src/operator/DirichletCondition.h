#ifndef DIRICHLETCONDITION_H
#define DIRICHLETCONDITION_H

#include "Operator.h"

/*!
 * \brief defines an identity operator for assembling a Dirichlet condition, the RHS should be set to the field one wishes to impose
 */

namespace hfox{

class DirichletCondition : public Operator{

  public:
    /*!
     * \brief constructor inheritance
     */
    using Operator::Operator;
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param detJacobians the values of the determinants of jacobian matrices at the integration points
     */
    void assemble(const std::vector< double > & detJacobians, const std::vector< EMatrix > & invJacobians);

};//DirichletCondition

}//hfox

#endif //DIRICHLETCONDITION_H
