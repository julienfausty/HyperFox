#ifndef FUNCTIONINTEGRATOR_H
#define FUNCTIONINTEGRATOR_H

#include <vector>
#include <string>
#include "Integrator.h"

namespace hfox{

class FunctionIntegrator : public Integrator{

  public:

    /*!
     * \brief use the parent constructor
     */
    using Integrator::Integrator;

    /*!
     * \brief useful typedef for defining the functions that come into play here
     */
    typedef std::function<void(std::vector<double>&, std::vector<double>*)> FieldFunc;

    /*!
     * \brief set the Function to integrate
     *
     * @params myField the function to actually integrate over
     */
    void setFunction(FieldFunc aFunc);

  protected:

    /*!
     * \brief evaluate the field values at the integration points of a given entity
     *
     * @param iEnt index of the entity to evaluate
     * @param ipVals array to fill with the values at the integration points
     */
    void evaluateIntegrand(int iEnt, std::vector<double> * ipVals);

    /*!
     * \brief the function to evaluate as the integrand over the mesh
     */
    FieldFunc myFunc;

};//FunctionIntegrator

};//hfox

#endif//FUNCTIONINTEGRATOR_H
