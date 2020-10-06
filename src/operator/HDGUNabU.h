#ifndef HDGUNABU_H
#define HDGUNABU_H

#include "HDGNonLinearOperator.h"
#include "Mass.h"

namespace hfox{

/*!
 * \brief Linearized operator for the u\nabla u term in the HDG method
 *
 * \f[
 * -\int_{e} (u^i_{(0)} u^j + u^i u^j_{(0)}) \nabla_i \varphi_j - \int_{e} (u^j_{0}\nabla_i u^{i} + u^j \nabla_i u^i_{(0)}) \varphi_j + \int_{\partial e} (u^i_{(0)} u^j + u^i u^j_{(0)}) n_i \varphi_j = -\int_{e} (u^i_{(0)} u^j_{(0)}) \nabla_i \varphi_j + \int_{\partial e} (u^i_{(0)} u^j_{(0)}) n_i \varphi_j
 * \f]
 */

class HDGUNabU : public HDGNonLinearOperator{

  public:
    /*!
     * \brief including default constructor
     *
     * @param re a pointer to the reference element
     */
    HDGUNabU(const ReferenceElement * re);
    /*!
     * \brief destructor for the class
     */
    ~HDGUNabU();
    /*!
     * \brief a method for allocating the operator
     *
     * @param nDOFsPerNode the number of degrees of freedom per node
     */
    void allocate(int nDOFsPerNodeUser);
    /*!
     * \brief method for assembling the operator
     *
     * @param invJacobians the values of the inverse jacobian matrices at the integration points
     * @param dV values of discrete mesure at the integration points
     */
    void assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians);
    /*!
     * \brief method for setting the solution vector
     *
     * @param sol a list of solution vectors per node
     */
    void setSolution(const std::vector<EVector> & sol);
    /*!
     * \brief a method that borrows objects from the HDGBase operator
     *
     * @param ns a pointer to the normals to the faces at their IPs
     */
    void setFromBase(const std::vector<EVector> * ns);
  protected:
    /*!
     * \brief the solution vector list at the IPs
     */
    std::vector<EVector> sols;
    /*!
     * \brief the solution at the face IPs
     */
    std::vector<EVector> faceSols;
    /*!
     * \brief the solution vector list at the nodes
     */
    std::vector<EVector> solNodes;
    /*!
     * \brief pointer to the list of normals to the faces at their IPs
     */
    const std::vector<EVector> * normals = NULL;

};//HDGUNabU

}//hfox

#endif //HDGUNABU_H
