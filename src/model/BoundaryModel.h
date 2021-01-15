#ifndef BOUNDARYMODEL_H
#define BOUNDARYMODEL_H

#include "FEModel.h"

namespace hfox{

/*!
 * \brief an enum with the different types of boundary models
 */
enum BoundaryModelType{
  /*!
   * \brief A boundary model implementing a condition uniquely on the solution
   */
  CGType,
  /*!
   * \brief A boundary model that implements a condition on the solution, flux and trace
   */
  HDGType
};

/*!
 * \brief An abstract class for models specific to boundaries
 */

class BoundaryModel : public FEModel{

  public:
    /*!
     * \brief constructor inheritance
     */
    using FEModel::FEModel;

    /*!
     * \brief simple get method that returns the BoundaryModelType
     */
    BoundaryModelType getBoundaryModelType() const {return myType;};

  protected:
    /*!
     * \brief define the type of boundary model
     */
    BoundaryModelType myType = CGType;

};//BoundaryModel

};//hfox

#endif//BOUNDARYMODEL_H
