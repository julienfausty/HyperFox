#ifndef HDGEMBEDDEDMODEL_H
#define HDGEMBEDDEDMODEL_H

#include "FEModel.h"
//#include "HDGEmbeddedBase.h"

namespace hfox{

/*!
 * \brief An interface class for HDG models meant to be used on lower dimensional meshes embedded in higher dimensional spaces
 */

class HDGEmbeddedModel : public FEModel{

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
    virtual void setFieldMap(const std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief allocating the local matrix and rhs
     */
    virtual void allocate(int nDOFsPerNode);
    /*!
     * \brief the HDG compute method
     */
    virtual void compute();
    /*!
     * \brief set the dimension of the embedding space
     */
    void setEmbeddingDimension(int dim);
  private:
    /*!
     * \brief initialize all the operators for the class
     */
    virtual void initializeOperators();
    /*!
     * \brief compute the local matrix
     */
    virtual void computeLocalMatrix()=0;
    /*!
     * \brief compute the local right hand side
     */
    virtual void computeLocalRHS()=0;
    /*!
     * \brief compute the geometrical information of the element (Jacobian matrices, metric tensors, Christoffel symbols)
     */
    void computeElementGeometry();
    /*!
     * \brief number of DOFs per node
     */
    int nDOFsPNode = 1;
    /*!
     * \brief the dimension of the embedding space
     */
    int embeddingDim = 0;
    /*!
     * \brief the values of the metric tensor at the integration points
     */
    std::vector<EMatrix> metricTensor;
    /*!
     * \brief the values of the inverse metric tensor at the integration points
     */
    std::vector<EMatrix> inverseMetricTensor;
    /*!
     * \brief the derivatives of metric tensor components (formated such that columns are metric tensor components and rows are derivatives in reference space)
     */
    std::vector<EMatrix> derivMetricTensor;
    /*!
     * \brief the computed Christoffel symbols at the integration points (christoffelSymbol[ip][contravariant][covariant, covariant])
     */
    std::vector< std::vector<EMatrix> > christoffelSymbols;
    /*!
     * \brief second derivatives of coordinates at the ips (hessians[ip][embeddingCoord][deriv, deriv])
     */
    std::vector< std::vector<EMatrix> > hessians;
    /*!
     * \brief normals to the element faces at the integration points
     */
    std::vector<EVector> normals;

};//HDGEmbeddedModel

};//hfox

#endif//HDGEMBEDDEDMODEL_H
