#ifndef UNABU_H
#define UNABU_H

#include "NonLinearOperator.h"
namespace hfox{

/*!
 * \brief Linearized operator for the u \nabla u term
 *
 * \f[
 * \int_{e} (u^i_{(0)} \nabla_i u^j + u^i \nabla_i u^j_{(0)}) \varphi_j = \int_{e} \nabla_i (u^i_{(0)} u^j_{(0)}) \varphi_j
 * \f]
 */

class UNabU : public NonLinearOperator{

  public:
    /*!
     * \brief constructor inheritance
     */
    using NonLinearOperator::NonLinearOperator;
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

};//UNabU

}//hfox

#endif //UNABU_H
