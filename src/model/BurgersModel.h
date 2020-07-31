#ifndef BURGERSMODEL_H
#define BURGERSMODEL_H

#include "FEModel.h"
#include "NabUU.h"
#include "Diffusion.h"

namespace hfox{

/*!
 * \brief A model implementing a linearized burgers equation in a CG setting
 *
 * \f[
 * \partial_{t} u^(i+1) + \nabla \cdot ( u^i \otimes u^(i+1) + u^(i+1) \otimes u^i - D \cdot \nabla u^(i+1))  = \nabla \cdot (u^i \otimes u^i)
 * \f]
 */

class BurgersModel : public FEModel{
  public:
    /*!
     * \brief constructor inheritance
     */
    using FEModel::FEModel;
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
};//BurgersModel

}//hfox

#endif//BURGERSMODEL_H
