#ifndef HDGCONVECTIONDIFFUSIONREACTIONSOURCE_H
#define HDGCONVECTIONDIFFUSIONREACTIONSOURCE_H

#include "HDGModel.h"
#include "HDGConvection.h"
#include "HDGDiffusion.h"
#include "Source.h"
#include "Reaction.h"

namespace hfox{

/*!
 * \brief A model implementing a convection-diffusion equation with possible reaction and source terms in an HDG setting
 *
 * \f[
 * \partial_{t} u + \nabla \cdot ( v u - D\nabla u) + r u  = f
 * \f]
 */

class HDGConvectionDiffusionReactionSource : public HDGModel{
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
    /*!
     * \brief set the source function
     */
    void setReactionFunction(std::function<double(const std::vector<double>&)> r);
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
     * \brief parse the velocity field values into vectors
     */
    std::vector<EVector> parseVelocityVals() const;
    /*!
     * \brief method for parsing diffusion tensor values into matrices
     */
    std::vector<EMatrix> parseDiffusionVals() const;
    /*!
     * \brief boolean tracking if a source function has been set
     */
    bool sourceSet = 0;
    /*!
     * \brief boolean tracking if reaction function has been set
     */
    bool reactionSet = 0;
};//HDGConvectionDiffusionReactionSource

}//hfox

#endif//HDGCONVECTIONDIFFUSIONREACTIONSOURCE
