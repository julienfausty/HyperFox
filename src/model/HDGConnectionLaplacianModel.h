#ifndef HDGCONNECTIONLAPLACIANMODEL_H
#define HDGCONNECTIONLAPLACIANMODEL_H

#include "HDGEmbeddedModel.h"
#include "Source.h"

namespace hfox{

/*!
 * \brief Implements a model for $-\nabla_{\Gamma} D \nabla_{\Gamma} u = \rho$ a connection laplacian Poisson problem.
 */

class HDGConnectionLaplacianModel : public HDGEmbeddedModel{

  public:

    /*!
     * \brief Constructor inheritance
     */
    using HDGEmbeddedModel::HDGEmbeddedModel;

    /*!
     * \brief a method to set the local fields
     *
     * @param pointer to a map of names and corresponding local values of fields
     */
    virtual void setFieldMap(const std::map<std::string, std::vector<double> > * fm);

    /*!
     * \brief the allocating of the matrix and the rhs
     */
    void allocate(int nDOFsPerNode);

    /*!
     * \brief set the source function
     *
     * @param s a function that associates a source to every coordinate point
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
    void parseDiffusionVals();
    /*!
     * \brief the values of the diffusion tensor at the integration points
     */
    std::vector<EMatrix> diffusionTensorVals;
    /*!
     * \brief boolean tracking if a source function has been set
     */
    bool sourceSet = 0;

};//HDGConnectionLaplacianModel

};//hfox

#endif//HDGCONNECTIONLAPLACIANMODEL_H
