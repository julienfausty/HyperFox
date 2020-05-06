#ifndef TIMESCHEME_H
#define TIMESCHEME_H

#include "Mass.h"

namespace hfox{

/*!
 * \brief TimeScheme is an interface to implementations of different time discretizations
 *
 * The idea here is to implement a time scheme as a speical mass operator with two supplemental methods:
 * - "setFieldMap": a method that takes the solution fields it needs from the provided field map
 * - "apply": a method that takes the stiffness matrix and rhs and can apply the time scheme to them
 */

class TimeScheme : public Mass{
  public:
    using Mass::Mass;
    /*!
     * \brief sets the field map for the class (needs to be virtual to specify certain name checks)
     *
     * @param fm the field map
     */
    virtual void setFieldMap(std::map<std::string, const std::vector<double> * > * fm)=0;
    /*!
     * \brief the method that applies the time scheme to the local stiffness matrix and rhs inplace
     *
     * @param stiffness pointer to the stiffness matrix that will become the total matrix
     * @param rhs pointer to the right hand side
     */
    virtual void apply(EMatrix * stiffness, EVector * rhs)=0;
    /*!
     * \brief method to set the time step
     *
     * @param timeStep the time step value
     */
    void setTimeStep(double timeStep){deltat = timeStep;};

  protected:
    /*!
     * \brief a pointer to the field map
     */
    std::map<std::string, const std::vector<double> * > * fieldMap = NULL;
    /*!
     * \brief time step parameter
     */
    double deltat = 0.0;
};

}//hfox

#endif//TIMESCHEME_H
