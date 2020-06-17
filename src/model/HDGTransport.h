#ifndef HDGTRANSPORT_H
#define HDGTRANSPORT_H

#include "HDGModel.h"
#include "HDGConvection.h"

namespace hfox{

/*!
 * \brief A model implementing the transport equation in an HDG setting
 *
 * \f[
 * \partial_{t} u - v \cdot \nabla u = 0
 * \f]
 */

class HDGTransport : public HDGModel{

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


};//HDGTransport

}//hfox


#endif //HDGTRANSPORT_H
