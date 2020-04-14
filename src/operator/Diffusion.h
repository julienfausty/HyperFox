#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <vector>
#include "Operator.h"
#include "ReferenceElement.h"

/*!
 * \brief The operator for diffusional terms
 *
 * \int_{Omega} \partial u \cdot \partial \varphi
 *
 */

namespace hfox{

class Diffusion : public Operator{

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

};//Diffusion

}//hfox

#endif //DIFFUSION_H
