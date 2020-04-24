#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <vector>
#include <algorithm>
#include "Operator.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief The operator for diffusional terms
 *
 * \int_{Omega} \partial u \cdot \partial \varphi
 *
 */

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
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the diffusion coefficient tensor
     *
     * @param diffCoeff a pointer to a list of diffusion tensors per node
     */
    void setDiffTensor(const std::vector<EMatrix> & diffCoeff);
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
