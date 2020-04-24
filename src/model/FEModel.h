#ifndef FEMODEL_H
#define FEMODEL_H

#include <map>
#include <string>
#include "Model.h"
#include "ReferenceElement.h"
#include "DenseEigen.h"
#include "ErrorHandle.h"
#include "Operator.h"
#include "AssemblyType.h"


namespace hfox{

/*!
 * \brief An interface class that describes the necessary elements of a finite element model
 *
 * The finite element model should implement the relevant physics towards the assembly of the linear system
 * for solving the finite element problem.
 */

class FEModel : public Model{
  public:
    /*!
     * \brief a method to set the reference element
     */
    FEModel(const ReferenceElement * re);
    /*!
     * \brief destructor
     */
    virtual ~FEModel();
    /*!
     * \brief allocate the local objects and initialize operators
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    virtual void allocate(int nDOFsPerNode) = 0;
    /*!
     * \brief the implementation of the compute method
     */
    virtual void compute();
    /*!
     * \brief a method to set the local fields
     *
     * @param pointer to a map of names and corresponding local values of fields
     */
    virtual void setFieldMap(const std::map<std::string, std::vector<double> > * fm) = 0;
    /*!
     * \brief compute the local matrix
     */
    virtual const EMatrix * getLocalMatrix() const{return &localMatrix;};
    /*!
     * \brief compute the local right hand side
     */
    virtual const EVector * getLocalRHS() const{return &localRHS;};
    /*!
     * \brief get the assembly type struct
     */
    virtual const AssemblyType * getAssemblyType() const{return &assembly;};
  protected:
    /*!
     * \brief a method where the relevant operators should be initialized (a good time to set the assembly type to)
     */
    virtual void initializeOperators()=0;
    /*!
     * \brief compute the local matrix
     */
    virtual void computeLocalMatrix()=0;
    /*!
     * \brief compute the local right hand side
     */
    virtual void computeLocalRHS()=0;
    /*!
     * \brief the reference element of the FE mesh
     */
    const ReferenceElement * refEl;
    /*!
     * \brief the local element matrix
     */
    EMatrix localMatrix;
    /*!
     * \brief the local element right hand side
     */
    EVector localRHS;
    /*!
     * \brief a map containing the local operators relevant to the model with their names
     */
    std::map<std::string, Operator*> operatorMap;
    /*!
     * \brief a vector of jacobians at IPs
     */
    std::vector<EMatrix> jacobians;
    /*!
     * \brief a vector of the inverse jacobians at IPs
     */
    std::vector<EMatrix> invJacobians;
    /*!
     * \brief a vector of the discrete measure values at IPs (detJac*weight)
     */
    std::vector<double> dV;
    /*!
     * \brief boolean value tracking if the local objects have been allocated
     */
    bool allocated = 0;
    /*!
     * \brief a structure defining if the matrix and rhs need to be set or added into the global objects
     */
    AssemblyType assembly;

}; //FEModel

}//hfox

#endif //FEMODEL_H
