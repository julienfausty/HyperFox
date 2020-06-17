#ifndef REACTION_H
#define REACTION_H

#include <vector>
#include <algorithm>
#include <functional>
#include "Mass.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief the operator implementing the reaction term from an analytical function
 *
 * \f[
 * \int_{\Omega} r u \varphi
 * \f]
 */

class Reaction : public Mass{
  public:
    /*!
     * \brief constructor inheritance
     */
    using Mass::Mass;
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the reaction function
     *
     * @param s the analytical reaction function (takes nodes)
     */
    void setReactionFunction(std::function<double(const std::vector<double>&)> r){reactionFunc = r;};
    /*!
     * \brief method for calculating source at IPs
     *
     * @param nodes element nodes
     */
    void calcReaction(const std::vector< std::vector<double> > & nodes);
  protected:
    /*!
     * \brief the source values at the IPs
     */
    std::vector<double> reaction;
    /*!
     * \brief the source function
     */
    std::function<double(const std::vector<double>&)> reactionFunc;
};//Reaction

}//hfox

#endif//REACTION_H
