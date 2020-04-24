#include <vector>
#include "Operator.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief The operator implementing the mass matrix of the element
 *
 * \int_{\Omega} u \varphi
 */

class Mass : public Operator{

  public:
    /*!
     * \brief constructor inheritance
     */
    using Operator::Operator;
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);

};//Mass

}//hfox
