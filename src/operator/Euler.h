#ifndef EULER_H
#define EULER_H

#include "TimeScheme.h"

namespace hfox{

/*!
 * \brief Euler class implemeting backward and forward Euler time integration strategies
 */

class Euler : public TimeScheme{

  public:
    /*!
     * \brief constructor
     */
    Euler(const ReferenceElement * re, bool isExplicitUser=false);
    /*!
     * \brief sets the field map for the class
     *
     * @param fm the field map
     */
    void setFieldMap(std::map<std::string, const std::vector<double> * > * fm);
    /*!
     * \brief the method that applies the time scheme to the local stiffness matrix and rhs inplace
     *
     * @param stiffness pointer to the stiffness matrix that will become the total matrix
     * @param rhs pointer to the right hand side
     */
    void apply(EMatrix * stiffness, EVector * rhs);
  protected:
    /*!
     * \brief boolean setting the behavior to explicit integration
     */
    bool isExplicit = false;

};//Euler

}//hfox

#endif//EULER_H
