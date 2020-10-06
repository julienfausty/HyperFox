#ifndef HDGBURGERSMODEL_H
#define HDGBURGERSMODEL_H

#include "HDGModel.h"
#include "HDGBase.h"
#include "HDGUNabU.h"
#include "HDGDiffusion.h"

namespace hfox{

/*!
 * \brief A model implementing a linearized burgers equation in an HDG setting
 *
 * \f[
 * \partial_{t} u^(i+1) + \nabla \cdot ( u^i \otimes u^(i+1) + u^(i+1) \otimes u^i - D \cdot \nabla u^(i+1))  = \nabla \cdot (u^i \otimes u^i)
 * \f]
 */

class HDGBurgersModel : public HDGModel{
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
     * \brief parse the solution field values into vectors
     */
    std::vector<EVector> parseSolutionVals() const;
    /*!
     * \brief method for parsing diffusion tensor values into matrices
     */
    std::vector<EMatrix> parseDiffusionVals() const;
};//HDGBurgersModel

}//hfox

#endif//HDGBURGERSMODEL_H
