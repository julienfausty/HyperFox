#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "FEModel.h"
#include "Convection.h"

namespace hfox{

/*!
 * \brief Continuous Galerkin model for the transport equation
 *
 * \[ \partial_t u + v \cdot \nabla u = 0 \]
 *
 */

class Transport : public FEModel{
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
     * \brief parse the velocity field values into vectors
     */
    std::vector<EVector> parseVelocityVals() const;
};//Transport

}//hfox

#endif//TRANSPORT_H
