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
    virtual void setBoundaryModel(FEModel * m){boundaryModel = m;};
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
     * \brief a helper method for constructing local fields
     *
     * @param entities indexes of the local entities
     * @param fm the field map to fill
     */
    void constructLocalFields(std::vector<int> & entities, std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief the interface to the linear algebra package
     */
    LinAlgebraInterface * linSystem = NULL;
    /*!
     * \brief the model of the system being simulated
     */
    FEModel * model = NULL;
    /*!
     * \brief the model of the boundary being simulated
     */
    FEModel * boundaryModel = NULL;
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
