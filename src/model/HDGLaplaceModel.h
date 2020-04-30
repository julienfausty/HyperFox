#ifndef HFGLAPLACEMODEL_H
#define HDGLAPLACEMODEL_H

#include <string>
#include <map>
#include <algorithm>
#include "FEModel.h"
#include "HDGBase.h"

namespace hfox{

/*!
 * \brief The model that implements the Laplace equation ($\Delta u = 0$) in an HDG setting
 */

class HDGLaplaceModel : public FEModel{

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
     * \brief compute the element jacobians, inverse jacobians and measures
     */
    void computeElementJacobians();
    /*!
     * \brief number of DOFs per node
     */
    int nDOFsPNode = 1;

};//HDGLaplaceModel

}

#endif//HDGLAPLACEMODEL_H
