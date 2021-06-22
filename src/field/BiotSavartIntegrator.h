#ifndef BIOTSAVARTINTEGRATOR_H
#define BIOTSAVARTINTEGRATOR_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <mpi.h>
#include "DenseEigen.h"
#include "Integrator.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief An Integrator class for performing Biot-Savart integrals
 */

class BiotSavartIntegrator : public Integrator{

  public:

    /*!
     * \brief use the parent constructor
     */
    BiotSavartIntegrator(Mesh * myMesh, int dim);

    /*!
     * \brief set the current field
     *
     * @params myCurrent a pointer to the current field on the mesh
     */
    void setCurrent(Field * myCurrent);

    /*!
     * \brief set the coordinates of the point to evaluate the magnetic field at
     *
     * @params coords a pointer to a vector of coordinates
     */
    void setCoordinates(std::vector<double> * coords);

    /*!
     * \brief set the dimension of the integral, always 3 for Biot Savart
     *
     * @param dim the dimension of the result of the integral, for Biot-Savart should always be equal to 3
     */
    void setDim(int dim);

    /*!
     * \brief main method of the class to perform the integral
     *
     * @param integral the return value in place (an array if it is multivalued)
     */
    void integrate(std::vector<double> * integral);

  protected:

    /*!
     * \brief evaluate the field values at the integration points of a given entity
     *
     * @param iEnt index of the entity to evaluate
     * @param ipVals array to fill with the values at the integration points
     */
    void evaluateIntegrand(int iEnt, std::vector<double> * ipVals);

    /*!
     * \brief current field
     */
    Field * pCurrent;

    /*!
     * \brief coordinates of the point to evaluate B at
     */
    std::vector<double> * pCoords;


};//BiotSavartIntegrator

};//hfox

#endif//BIOTSAVARTINTEGRATOR_H
