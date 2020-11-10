#ifndef HDGNGAMMAMODEL_H
#define HDGNGAMMAMODEL_H

#include "nGammaParams.h"
#include "HDGModel.h"
#include "HDGBase.h"
#include "HDGDiffusion.h"
#include "Source.h"

namespace hfox{

namespace nGamma{

/*!
 * \brief A model implementing a linearized nGamma model
 *
 * \f[
 *    \left\{
 *    \begin{array}{l}
 *    \partial_t n + \nabla_{\perp}(\Gamma b + nu_{\not b} - D\nabla_{\perp}n) = S_{0}\\
 *    \partial_t (\Gamma) + \nabla_{\perp} \cdot \left(\dfrac{\Gamma^{2}}{n}b + \Gamma u_{\not b} + nc_{s}^{2}b - G\nabla_{\perp}\Gamma\right) = S_1
 *    \end{array}
 *    \right .
 * \f]
 */

class HDGnGammaModel : public HDGModel {
  public:
    /*!
     * \brief constructor inheritance
     */
    using HDGModel::HDGModel;
    /*!
     * \brief a method to set the local fields
     *
     * @param pointer to a map of names and corresponding local values of fields
     */
    void setFieldMap(const std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief allocating the local matrix and rhs
     */
    void allocate(int nDOFsPerNode);
    /*!
     * \brief set the source function
     */
    void setSourceFunction(std::function<double(const std::vector<double>&, int)> s);
  protected:
    /*!
     * \brief initialize all the operators for the class
     */
    void initializeOperators();
    /*!
     * \brief compute the local matrix
     */
    void computeLocalMatrix();
    /*!
     * \brief compute the local right hand side
     */
    void computeLocalRHS();
    /*!
     * \brief specific parameter struct for nGamma models
     */
    nGammaParams params;
    /*!
     * \brief boolean that tracks if the source has been set
     */
    bool sourceSet = 0;
};//HDGnGammaModel

};//nGamma

}//hfox

#endif//HDGNGAMMAMODEL_H
