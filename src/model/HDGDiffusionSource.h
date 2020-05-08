#ifndef HDGDIFFUSIONSOURCE_H
#define HDGDIFFUSIONSOURCE_H

#include "HDGModel.h"
#include "HDGBase.h"
#include "HDGDiffusion.h"
#include "TimeScheme.h"
#include "Source.h"

namespace hfox{

/*!
 * \brief A model implementing the diffusion with a source term equation in an HDG setting
 *
 * \[
 * \partial_{t} u - \Delta u = f
 * \]
 */

class HDGDiffusionSource : public HDGModel{

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

};//HDGDiffusionSource

}//hfox

#endif//HDGDIFFUSIONSOURCE_H
