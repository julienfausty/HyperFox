#ifndef HDGBASE_H
#define HDGBASE_H

#include "ErrorHandle.h"
#include "DenseEigen.h"
#include "HDGOperator.h"
#include "Mass.h"
#include "Convection.h"

namespace hfox{

/*!
 * \brief base operator for the HDG method (the terms that are always present)
 *
 * \f[
 * \int_{\partial e} -\tau(u-\lambda)\varphi
 * \f]
 * \int_{e} q\omega + \int_{e} u\div(\omega) - \int_{\partial e}u (\omega\cdot n)
 * \f[
 * \f]
 * \int_{\partial e} (q \cdot n) \eta + \int_{\partial e} \tau(u-\lambda)\eta
 * \f[
 * \f]
 */

class HDGBase : public HDGOperator{

  public:
    /*!
     * \brief including default constructors
     */
    HDGBase(const ReferenceElement * re);
    /*!
     * \brief destructor for the class
     */
    ~HDGBase();
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
     * \brief method for setting the stabilization field tau
     *
     * @param tauUser a list of stabilization terms per node
     */
    void setTau(const std::vector<double> & userTaus);
    /*!
     * \brief method for calculating normals at the IPs
     *
     * @param the nodes of the element
     * @param the jacobians at the IPs of the element and all the face elements
     */
    void calcNormals(const std::vector< std::vector<double> > & nodes, const std::vector<EMatrix> & jacobians);
    /*!
     * \brief get method for bulk mass operator
     */
    const Mass * getBulkMass() const{return bulkMass;};
    /*!
     * \brief get method for face mass operator
     */
    const Mass * getFaceMass() const{return faceMass;};
    /*!
     * \brief get method for face mass operator
     */
    const Convection * getConvection() const{return convection;};
  protected:
    /*!
     * \brief the stabilization term list at the IPs
     */
    std::vector<double> taus;
    /*!
     * \brief the normals list at the IPs
     */
    std::vector<EVector> normals;
    /*!
     * \brief a face mass operator
     */
    Mass * faceMass;
    /*!
     * \brief a bulk mass operator
     */
    Mass * bulkMass;
    /*!
     * \brief a bulk convection operator
     */
    Convection * convection;

};//HDGBase

}//hfox

#endif//HDGBASE_H
