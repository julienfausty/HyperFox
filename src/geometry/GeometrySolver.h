#ifndef GEOMETRYSOLVER_H
#define GEOMETRYSOLVER_H

#include <vector>
#include "Mesh.h"
#include "Partitioner.h"
#include "DenseEigen.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief An interface for solvers that compute aspects of the mesh geometry
 *
 * Using the finite element description of the mesh, this type of solver should calculate any geometrical information related to it (i.e. metric
 * tensor, parameterized coordinates, jacobians, etc.)
 */

class GeometrySolver {

  public:
    /*!
     * \brief The default constructor that just takes a mesh
     *
     * @param pmesh a pointer to a mesh object
     */
    GeometrySolver(const Mesh * pmesh){myMesh = pmesh;};

    //Destructor
    virtual ~GeometrySolver(){};

    /*!
     * \brief Virtual method for running the geometry solver
     */
    virtual void solve() = 0;

  protected:

    const Mesh * myMesh;

};//GeometrySolver

};//hfox

#endif//GEOMETRYSOLVER_H
