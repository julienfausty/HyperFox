#ifndef HDGMODEL_H
#define HDGMODEL_H

#include <string>
#include <map>
#include <algorithm>
#include "FEModel.h"
#include "HDGBase.h"

namespace hfox{

/*!
 * \brief Interface for HDG models
 */

class HDGModel : public FEModel{
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
  protected:
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
     * \brief compute the element jacobians, inverse jacobians and measures
     */
    void computeElementJacobians();
    /*!
     * \brief number of DOFs per node
     */
    int nDOFsPNode = 1;

};//HDGModel

}//hfox

#endif//HDGMODEL_H
