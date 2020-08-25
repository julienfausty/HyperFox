#ifndef HDGNABUU_H
#define HDGNABUU_H

#include "HDGNonLinearOperator.h"

namespace hfox{

/*!
 * \brief Linearized operator for the \nabla (u \otimes u) term in the HDG method
 *
 * \f[
 * -\int_{e} (u^i_{(0)} u^j + u^i u^j_{(0)}) \nabla_i \varphi_j + \int_{\partial e} (u^i_{(0)} u^j + u^i u^j_{(0)}) n_i \varphi_j = -\int_{e} (u^i_{(0)} u^j_{(0)}) \nabla_i \varphi_j + \int_{\partial e} (u^i_{(0)} u^j_{(0)}) n_i \varphi_j
 * \f]
 */

class HDGNabUU : public HDGNonLinearOperator{

  public:
    /*!
     * \brief including default constructor
     *
     * @param re a pointer to the reference element
     */
    HDGNabUU(const ReferenceElement * re);
    /*!
     * \brief destructor for the class
     */
    ~HDGNabUU();
    /*!
     * \brief a method for allocating the operator
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    void allocate(int nDOFsPerNodeUser);
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the solution vector
     *
     * @param sol a list of solution vectors per node
     */
    void setSolution(const std::vector<EVector> & sol);
  protected:
    /*!
     * \brief the solution vector list at the IPs
     */
    std::vector<EVector> sols;
    /*!
     * \brief the solution vector list at the nodes
     */
    std::vector<EVector> solNodes;
    /*!
     * \brief the gradient of the solution at the IPs
     */
    std::vector<EMatrix> gradSols;

};//HDGNabUU

}//hfox

#endif //HDGNABUU_H
