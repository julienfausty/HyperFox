#ifndef OPERATOR_H
#define OPERATOR_H

#include "Field.h"
#include "DenseEigen.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief an object representing an operator on a field
 *
 * This object should take a vector of values defined in an element to another vector of values defined in the same element.
 */

class Operator{
  public:
    //Default constructor
    Operator(ReferenceElement * referenceEl){refEl = referenceEl;};
    // Destructor
    virtual ~Operator(){};
    /*!
     * \brief access the operator matrix
     */
    const EMatrix * getMatrix() const {return &op;};
    /*!
     * \brief a virtual method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param detJacobians the values of the determinants of jacobian matrices at the integration points
     */
    virtual void assemble(const std::vector< double > & detJacobians, 
        const std::vector< EMatrix > & invJacobians) = 0;
    /*!
     * \brief a method for allocating the operator
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    virtual void allocate(int nDOFsPerNodeUser);
    /*!
     * \brief a helper method for computing the Jacobian at the integration points of an element
     *
     * @param points nodes of the mesh element
     * @param referenceEl a pointer to the reference element
     */
    static std::vector<EMatrix> calcJacobians(const std::vector< std::vector<double> > & points, 
        const ReferenceElement * referenceEl);
  protected:
    /*!
     * \brief the matrix representation of the operator
     */
    EMatrix op;
    /*!
     * \brief the reference element
     */
    const ReferenceElement * refEl;
    /*!
     * \brief the number of degrees of freedom per node
     */
    int nDOFsPerNode;
    /*!
     * \brief a boolean tracking the allocation of the operator
     */
    bool allocated = 0;
};//Operator

} //hfox

#endif
