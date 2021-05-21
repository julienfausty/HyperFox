#ifndef FIELDINTEGRATOR_H
#define FIELDINTEGRATOR_H

#include <vector>
#include <mpi.h>
#include "DenseEigen.h"
#include "Integrator.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief An Integrator class for integrating over fields
 */

class FieldIntegrator : public Integrator{

  public:

    /*!
     * \brief use the parent constructor
     */
    using Integrator::Integrator;

    /*!
     * \brief set the Field to integrate
     *
     * @params myField a pointer to the field to actually integrate over
     */
    void setField(Field * myField);

  protected:

    /*!
     * \brief evaluate the field values at the integration points of a given entity
     *
     * @param iEnt index of the entity to evaluate
     * @param ipVals array to fill with the values at the integration points
     */
    void evaluateIntegrand(int iEnt, std::vector<double> * ipVals);

    /*!
     * \brief field to integrate
     */
    Field * pField;


};//FieldIntegrator

};//hfox

#endif//FIELDINTEGRATOR_H
