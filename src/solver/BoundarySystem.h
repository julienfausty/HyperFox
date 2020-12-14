#ifndef BOUNDARYSYSTEM_H
#define BOUNDARYSYSTEM_H

#include <set>
#include <vector>
#include <tuple>
#include "FEModel.h"

namespace hfox{

/*!
 * \brief A boundary condition manager
 *
 * BoundarySystem is a manager class for boundary conditions. It is meant to do the boundary book keeping 
 * for the Solver class.
 */

typedef std::vector< std::tuple<FEModel*, const std::set<int> * > > BoundaryConditionList;

class BoundarySystem{

  public:

    /*!
     * \brief empty constructor
     */
    BoundarySystem(){};

    /*!
     * \brief method to set a boundary condition
     *
     * @param bModel a model to set on the boundary as a boundary condition
     * @param faces a pointer to a set of face indeces that constitute the boundary
     */
    void setBoundaryCondition(FEModel * bModel, const std::set<int> * faces);

    /*!
     * \brief a method to get a constant pointer to the boundary list
     */
    const BoundaryConditionList * getBoundaryList(){return &boundaryList;};

  private:

    /*!
     * \brief data structure for holding boundary model/face set pairs
     */
    BoundaryConditionList boundaryList;

};//BoundarySystem

};//hfox

#endif///BOUNDARYSYSTEM_H
