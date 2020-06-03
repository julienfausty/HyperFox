#ifndef HFGLAPLACEMODEL_H
#define HDGLAPLACEMODEL_H

#include <string>
#include <map>
#include <algorithm>
#include "HDGModel.h"
#include "HDGBase.h"
#include "HDGDiffusion.h"

namespace hfox{

/*!
 * \brief The model that implements the Laplace equation ($\Delta u = 0$) in an HDG setting
 */

class HDGLaplaceModel : public HDGModel{

  public:
    /*!
     * \brief constructor inheritance
     */
    using HDGModel::HDGModel;
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

};//HDGLaplaceModel

}

#endif//HDGLAPLACEMODEL_H
