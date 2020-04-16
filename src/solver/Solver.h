#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include "LinAlgebraInterface.h"
#include "FEModel.h"
#include "Field.h"
#include "Mesh.h"

namespace hfox{

/*!
 * \brief The interface to create solvers
 *
 * The main goal of this interface to to be able to assemble and solve systems 
 * independently of the models that are being used.
 */

class Solver{
  public:
    /*!
     * \brief a method to set the model
     */
    virtual void setModel(FEModel * m){model = m;};
    /*!
     * \brief a method to set the model
     */
    virtual void setBoundaryModel(FEModel * m){boundaryModel = m;};
    /*!
     * \brief a method to set the linear algebra interface
     */
    virtual void setLinSystem(LinAlgebraInterface * lai){LinSystem = lai;};
    /*!
     * \brief a method to set the linear algebra interface
     */
    virtual void setFieldMap(std::map<std::string, Field* > * fm){fieldMap = fm;};
    /*!
     * \brief a method to set the linear algebra interface
     */
    virtual void setMesh(Mesh * m){myMesh = m;};
    /*!
     * \brief a method to assemble the linear system
     */
    virtual void assemble()=0;
    /*!
     * \brief a method to solve the linear system
     */
    virtual void solve()=0;
  protected:
    /*!
     * \brief the interface to the linear algebra package
     */
    LinAlgebraInterface * LinSystem;
    /*!
     * \brief the model of the system being simulated
     */
    FEModel * model;
    /*!
     * \brief the model of the boundary being simulated
     */
    FEModel * boundaryModel;
    /*!
     * \brief the map to the fields
     */
    std::map<std::string, Field* > * fieldMap;
    /*!
     * \brief a pointer to the mesh
     */
    Mesh * myMesh;

};//Solver

}//hfox

#endif //SOLVER_H
