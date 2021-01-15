#ifndef INTEGRATEDDIRICHLETMODEL_H
#define INTEGRATEDDIRICHLETMODEL_H

#include "BoundaryModel.h"
#include "Mass.h"

namespace hfox{

/*!
 * \brief A model that implements Dirichlet boundary conditions: 
 * \[
 * \int u \varphi = \int u_d \varphi 
 * \]
 * in Element
 */

class IntegratedDirichletModel : public BoundaryModel{

  public:
    /*!
     * \brief constructor inheritance
     */
    using BoundaryModel::BoundaryModel;
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

};//IntegratedDirichletModel

}//hfox

#endif//INTEGRATEDDIRICHLETMODEL_H
