#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include "LinAlgebraInterface.h"
#include "FEModel.h"
#include "BoundaryModel.h"
#include "BoundarySystem.h"
#include "Field.h"
#include "Mesh.h"
#include "Partitioner.h"
#include "Utils.h"

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
     * \brief empty constructor
     */
    Solver(){};
    /*!
     * \brief a method to set the model
     *
     * @param m a pointer to the model
     */
    virtual void setModel(FEModel * m){model = m;};
    /*!
     * \brief a method to set the boundary model
     *
     * @param m a pointer to the boundary model
     */
    virtual void setBoundaryModel(BoundaryModel * m){bSystem.setBoundaryCondition(m, myMesh->getBoundaryFaces());};
    /*!
     * \brief a method to set the boundary condition on a specified set of faces
     *
     * @param m a pointer to the boundary model
     * @param faces a pointer to a set of faces
     */
    virtual void setBoundaryCondition(BoundaryModel * m, std::set<int> * faces){bSystem.setBoundaryCondition(m, faces);};
    /*!
     * \brief a method to set the linear algebra interface
     *
     * @param lai a pointer to the linear algebra interface
     */
    virtual void setLinSystem(LinAlgebraInterface * lai){linSystem = lai;};
    /*!
     * \brief a method to set the field map
     *
     * @param fm pointer to the field map
     */
    virtual void setFieldMap(std::map<std::string, Field* > * fm){fieldMap = fm;};
    /*!
     * \brief a method to set the mesh
     *
     * @param m pointer to the mesh
     */
    virtual void setMesh(Mesh * m){myMesh = m;};
    /*!
     * \brief a method to set the verbosity (0 quiet | 1 verbose)
     *
     * @param v boolean controlling verbosity
     */
    virtual void setVerbosity(bool v){verbose = v;};
    /*!
     * \brief a method for initializing all components
     */
    virtual void initialize();
    /*!
     * \brief a method for allocating the memory for all components
     */
    virtual void allocate()=0;
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
     * \brief calculate the sparsity pattern of the matrix based on the mesh connectivity
     */
    virtual void calcSparsityPattern()=0;
    /*!
     * \brief a helper method to create a local field map with a given field type
     *
     * @param ft the field type the fieldMap is adhereing to
     */
    std::map<std::string, std::vector<double> > prepareLocalFieldMap(FieldType ft);
    /*!
     * \brief a helper method for constructing local fields (sequential)
     *
     * @param entities indexes of the local entities
     * @param fm the field map to fill
     */
    void constructLocalFields(std::vector<int> & entities, std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief a helper method for constructing local fields (parallel)
     *
     * @param entities global indexes of the entities
     * @param locentities local indexes of the entities
     * @param fm the field map to fill
     */
    void constructLocalFields(std::vector<int> & entities, std::vector<int> locEntities, std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief a helper method for reordering a vector in parallel
     * @param solutionVec pointer to the vector
     * @param thisRange the current indexes of the DOFs of the solution
     * @param dofsWeNeed the indexes this partition requires
     */
    virtual void reorder(std::vector<double> * solutionVec, std::vector<int> * thisRange, std::vector<int> * dofsWeNeed);
    /*!
     * \brief the interface to the linear algebra package
     */
    LinAlgebraInterface * linSystem = NULL;
    /*!
     * \brief the model of the system being simulated
     */
    FEModel * model = NULL;
    /*!
     * \brief the boundary system being simulated
     */
    BoundarySystem bSystem;
    /*!
     * \brief the map to the fields
     */
    std::map<std::string, Field* > * fieldMap;
    /*!
     * \brief a pointer to the mesh
     */
    Mesh * myMesh = NULL;
    /*!
     * \brief the diagonal number of non zero entries in the fe matrix per row
     */
    std::vector<int> diagSparsePattern;
    /*!
     * \brief the off-diagonal number of non zero entries in the fe matrix per row
     */
    std::vector<int> offSparsePattern;
    /*!
     * \brief the number of degrees of freedom per node
     */
    int nDOFsPerNode = 1;
    /*!
     * \brief boolean tracking initialization
     */
    bool initialized = 0;
    /*!
     * \brief boolean tracking allocation
     */
    bool allocated = 0;
    /*!
     * \brief boolean tracking assembly
     */
    bool assembled = 0;
    /*!
     * \brief boolean value controlling verbosity
     */
    bool verbose = 1;

};//Solver

}//hfox

#endif //SOLVER_H
