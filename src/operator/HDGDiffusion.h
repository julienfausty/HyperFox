#ifndef HDGDIFFUSION_H
#define HDGDIFFUSION_H

#include "HDGOperator.h"
#include "Convection.h"
#include "Mass.h"

namespace hfox{

/*!
 * \brief operator that implements diffusion in an HDG setting
 *
 * \[f
 * \int_{e} q \cdot \nabla \varphi - \int_{\partial e} q \cdot n \varphi
 * \]f
 *
 * since this operator uses a lot of the same components as HDGBase, instead of re allocating and 
 * recalculating them, this operator takes them a set values from the method "setFromBase"
 */

class HDGDiffusion : public HDGOperator{

  public:
    /*!
     * \brief including default constructor
     */
    HDGDiffusion(const ReferenceElement * re);
    /*!
     * \brief destructor for the class
     */
    ~HDGDiffusion();
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
     * \brief method for setting the diffusion tensor field
     *
     * @param diffTensor a list of diffusion tensors at the nodes of the element
     */
    void setDiffusionTensor(const std::vector<EMatrix> & diffTensor);
    /*!
     * \brief a method that borrows objects from the HDGBase operator
     *
     * @param cv a pointer to a pre allocated convection operator
     * @param fcmass a pointer to a pre allocated face mass operator
     * @param a pointer to the normals to the faces at their IPs
     */
    void setFromBase(const std::vector<EVector> * ns);

  protected:
    /*!
     * \brief list of diffusion tensors at the face IPs
     */
    std::vector<EMatrix> Ds;
    /*!
     * \brief list of diffusion tensors at the element nodes
     */
    std::vector<EMatrix> myDiffTensor;
    /*!
     * \brief a pointer to a convection operator
     */
    Convection * convection = NULL;
    /*!
     * \brief a pointer to a face mass operator
     */
    Mass * faceMass = NULL;
    /*!
     * \brief pointer to the list of normals to the faces at their IPs
     */
    const std::vector<EVector> * normals = NULL;
};//HDGOperator

}//hfox

#endif//HDGDIFFUSION_H
