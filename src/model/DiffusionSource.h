#ifndef DIFFUSIONSOURCE_H
#define DIFFUSIONSOURCE_H

#include "FEModel.h"
#include "Diffusion.h"
#include "Source.h"

namespace hfox{

/*!
 * \brief Continuous Galerkin model for the diffusion source equation
 *
 * \[
 * \partial_t u - D \Delta u = f
 * \]
 */

class DiffusionSource : public FEModel{
  
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
    /*!
     * \brief set the source function
     */
    void setSourceFunction(std::function<double(const std::vector<double>&)> s);
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
     * \brief method for parsing diffusion tensor values into matrices
     */
    std::vector<EMatrix> parseDiffusionVals() const;

};//DiffusionSource

};//hfox

#endif//DIFFUSIONSOURCE_H
