#ifndef SOURCE_H
#define SOURCE_H

#include <functional>
#include "RHSOperator.h"

namespace hfox{

/*!
 * \brief the operator implementing the source term from an analytical function
 *
 * \[\
 * \int_{\Omega} f \varphi
 * ]
 */

class Source : public RHSOperator{

  public:
    /*!
     * \brief constructor inheritance
     */
    using RHSOperator::RHSOperator;
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the source function
     *
     * @param s the analytical source function (takes nodes)
     */
    void setSourceFunction(std::function<double(const std::vector<double>&)> s){sourceFunc = s;};
    /*!
     * \brief method for calculating source at IPs
     *
     * @param nodes element nodes
     */
    void calcSource(const std::vector< std::vector<double> > & nodes);
  protected:
    /*!
     * \brief the source values at the IPs
     */
    std::vector<double> source;
    /*!
     * \brief the source function
     */
    std::function<double(const std::vector<double>&)> sourceFunc;

};//Source

}//hfox

#endif //SOURCE_H
