#ifndef CONVECTION_H
#define CONVECTION_H

#include "Operator.h"
#include "DenseEigen.h"
#include "ErrorHandle.h"

namespace hfox{

/*!
 * \brief The operator for convection terms
 *
 * \f[
 * \int_{\Omega} v \cdot \partial u \varphi 
 * \f]
 */

class Convection : public Operator{

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
     * \brief method for setting the velocity vector
     *
     * @param velocity a list of velocity vectors per node
     */
    void setVelocity(const std::vector<EVector> & velocity);
  protected:
    /*!
     * \brief the velocity vector list at the IPs
     */
    std::vector<EVector> vels;

};//Convection

}//hfox

#endif//CONVECTION_H
