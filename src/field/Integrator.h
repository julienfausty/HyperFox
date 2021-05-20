#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <mpi.h>
#include "Mesh.h"
#include "Field.h"
#include "Partitioner.h"

namespace hfox{

/*!
 * \brief Interface class for computing integrals on meshes.
 */

class Integrator{

  public:

    /*!
     * \brief main method of the class to perform the integral
     *
     * @param integral the return value in place (an array if it is multivalued)
     */
    virtual void integrate(std::vector<double> * integral) =0;

};//Integrator

};//hfox

#endif//INTEGRATOR_H
