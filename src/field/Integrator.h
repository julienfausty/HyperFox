#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <mpi.h>
#include "DenseEigen.h"
#include "Mesh.h"
#include "Partitioner.h"
#include "Operator.h"

namespace hfox{

/*!
 * \brief Interface class for computing integrals on meshes.
 */

class Integrator{

  public:

    /*!
     * \brief simplest constructor of the class
     *
     * @param myMesh pointer to a mesh
     * @param dim the dimension of the integral
     */
    Integrator(Mesh * myMesh, int dim){
      setMesh(myMesh);
      setDim(dim);
    };

    /*!
     * \brief main method of the class to perform the integral
     *
     * @param integral the return value in place (an array if it is multivalued)
     */
    void integrate(std::vector<double> * integral);

    /*!
     * \brief a set method for seting the mesh
     *
     * @param myMesh a pointer towards a Mesh object
     */
    void setMesh(Mesh * myMesh){
      pmesh = myMesh;
      if(pmesh == NULL){
        throw(ErrorHandle("FieldIntegrator", "setMesh", "cannot set a mesh as NULL."));
      }
    };

    /*!
     * \brief a get method to look at the mesh one is integrating over
     */
    Mesh * getMesh() const{return pmesh;};

    /*!
     * \brief a method for setting the dimension of the integral
     *
     * @param dim the dimension of the outgoing integral
     */
    void setDim(int dim){
      dimIntegral = dim;
      if(dimIntegral <= 0){
        throw(ErrorHandle("FieldIntegrator", "setDim", "cannot set the integral dimension to 0 or under."));
      }
    };

    /*!
     * \brief a get method for the dimension of the integral
     */
    int getDim() const{return dimIntegral;};

    /*!
     * \brief a set method for seting the type of entity to integrate over
     *
     * @param myType the type of entity to integrate over
     */
    void setType(FieldType myType){
      type = myType;
      if(type != Cell and type != Face){
        throw(ErrorHandle("Integrator", "setType", "only allowed to set entity types of Cell or Face."));
      }
    };

  protected:

    /*!
     * \brief the method to evaluate the integrand at the integration points
     *
     * @param iEnt the index of the entity one is integrating over
     * @param ipVals an array with the values of the integrand at the integration points
     */
    virtual void evaluateIntegrand(int iEnt, std::vector<double> * ipVals) =0;

    /*!
     * \brief the mesh one is integrating over
     */
    Mesh * pmesh = NULL;

    /*!
     * \brief the number of values in the integral
     */
    int dimIntegral = 0;

    /*!
     * \brief the entity type to perform the integral over (Face or Cell for now)
     */
    FieldType type = Cell;

};//Integrator

};//hfox

#endif//INTEGRATOR_H
