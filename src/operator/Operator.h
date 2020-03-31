#ifndef OPERATOR_H
#define OPERATOR_H

#include "Field.h"
#include "DenseEigen.h"
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief an object representing an operator on a field
 *
 * This object should take a field to another field in the same space as the original field.
 */

class Operator{
  public:
    // Destructor
    virtual ~Operator();
    /*!
     * \brief access the operator matrix
     */
    const EMatrix * getMatrix() const;
    /*!
     * \brief a virtual method for assembling the operator
     *
     * @param nodes are the values of the nodes of the element.
     */
    virtual void assemble(const std::vector< EVector > * nodes) = 0;
  protected:
    /*!
     * \brief a helper method for computing the Jacobian at a point in an element
     *
     * @param point the point in the element where the Jacobian should be calculated
     */
    EMatrix calcJacobian(const EVector & point) const;
    /*!
     * \brief the matrix representation of the operator
     */
    EMatrix op;
    /*!
     * \brief the reference element
     */
    const ReferenceElement * refEl;
};//Operator

} //hfox

#endif
