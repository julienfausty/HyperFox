#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <vector>
#include <algorithm>
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
    /*!
     * \brief method for setting the diffusion coefficient tensor
     *
     * @param diffCoeff a pointer to a list of diffusion tensors per node
     */
    void setDiffTensor(std::vector<EMatrix> & diffCoeff);
  protected:
    /*!
     * \brief helper method for extracting n(n+1)/2 vector from a symmetric matrix
     *
     * @param symMat a symmetric matrix
     */
    EVector symMat2Vec(const EMatrix & symMat);
    /*!
     * \brief a vector with the diffusion coefficient tensors at the IPs
     */
    std::vector<EMatrix> diffTensor;

};//Diffusion

}//hfox

#endif //DIFFUSION_H
