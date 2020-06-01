#ifndef HDGCONVECTION_H
#define HDGCONVECTION_H

#include "HDGOperator.h"
#include "Mass.h"
#include "Convection.h"

namespace hfox{

/*!
 *  \brief Convection operator for hybrid discontinuous Galerkin solver
 *
 *  \[
 *  \int_{e} v \cdot q \varphi
 *  \]
 */

class HDGConvection : public HDGOperator{
  public:
    /*!
     * \brief including default constructor
     *
     * @param re a pointer to the reference element
     */
    HDGConvection(const ReferenceElement * re);
    /*!
     * \brief destructor for the class
     */
    ~HDGConvection();
    /*!
     * \brief a method for allocating the operator
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    void allocate(int nDOFsPerNodeUser);
    /*!
     * \brief a virtual method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of the discrete measure
     */
    virtual void assemble(const std::vector< double > & dV, 
        const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the velocity field
     *
     * @param vels a list of velocities at the nodes of the element
     */
    void setVelocity(const std::vector<EVector> & vels);
  protected:
    /*!
     * \brief the velocities at the IPs
     */
    std::vector<EVector> velocities;
    /*!
     * \brief a pointer to the mass operator
     */
    Mass * mass = NULL;

};//HDGConvection

}//hfox

#endif//HDGCONVECTION_H
