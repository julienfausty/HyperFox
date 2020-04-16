#ifndef DIRICHLETMODEL_H
#define DIRICHLETMODEL_H

#include "FEModel.h"

namespace hfox{

/*!
 * \brief A model that implements Dirichlet boundary conditions: u = u_d in Element
 */

class DirichletModel : public FEModel{

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

};//DirichletModel

}//hfox

#endif//DIRICHLETMODEL_H
