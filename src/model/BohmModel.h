#ifndef BOHMMODEL_H
#define BOHMMODEL_H

#include <functional>
#include "FEModel.h"
#include "Mass.h"

namespace hfox{

namespace nGamma{

/*!
 * \brief A model that implements a Bohm-Chodura like boundary for an nGamma system
 */

class BohmModel : public FEModel{

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
     * \brief set the transfer function
     */
    void setTransferFunction(std::function<double(double, double)> tf, std::function<std::vector<double>(double, double)> derivTf){
      transferFunction = tf; 
      derivTransferFunction = derivTf;
    };
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
     * \brief transfer function
     */
    std::function<double(double, double)> transferFunction;
    /*!
     * \brief the vector derivative of the transfer function
     */
    std::function<std::vector<double>(double, double)> derivTransferFunction;
    /*!
     * \brief originalNonDerived part of the system
     */
    EMatrix originalSystem;

};//BohmModel

};//nGamma

};//hfox

#endif//BOHMMODEL_H
