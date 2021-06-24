#ifndef IP2NODESOLVER_H
#define IP2NODESOLVER_H

#include <functional>
#include <vector>
#include "GeometrySolver.h"
#include "Operator.h"
#include "Mass.h"
#include "Field.h"
#include "DenseEigen.h"
#include "ErrorHandle.h"

namespace hfox{

/*!
 * \brief A geometry solver that takes a field defined on the integration points and solves for it on the nodes of the mesh
 */

class IP2NodeSolver : public GeometrySolver{

  public:
    /*!
     * \brief Constructor inheritance
     */
    using GeometrySolver::GeometrySolver;

    /*!
     * \brief set the field to store the values on the nodes (! must be a Cell field)
     */
    void setField(Field * pField);

    /*!
     * \brief set function for the computation of values on integration points given the coordinates of the said integration points
     */
    void setFunction(std::function<void(std::vector<double>&,std::vector<double>*)> ipFunc);

    /*!
     * \brief set function for the computation of values on integration points given the jacobians at the said integration points
     */
    void setFunction(std::function<void(std::vector<EMatrix>&,std::vector<double>*)> ipFunc);

    /*!
     * \brief set function for the computation of values on integration points given the coordinates index of the cell
     */
    void setFunction(std::function<void(int,std::vector<double>*)> ipFunc);

    /*!
     * \brief the solve operation (F = M^{-1}b)
     */
    void solve();

    //Utility functions to set as functions of the class for easy computation
    /*!
     * \brief for returning all the jacobian matrix values in a single vector
     */
    static void jacobianVals(std::vector<EMatrix> & jacs, std::vector<double> * vals);
    /*!
     * \brief for returning computing the Riemannian metric values of the mesh
     */
    static void metricVals(std::vector<EMatrix> & jacs, std::vector<double> * vals);
  protected:

    /*!
     * \brief the field to put the values in
     */
    Field * myField;

    /*!
     * \brief the function from nodes to values
     */
    std::function< void(std::vector<double>&, std::vector<double>*) > ipFuncNodes;

    /*!
     * \brief the function from jacobians to values
     */
    std::function< void(std::vector<EMatrix>&, std::vector<double>*) > ipFuncJacs;

    /*!
     * \brief the function from cell index to values
     */
    std::function< void(int, std::vector<double>*) > ipFuncCell;

    /*!
     * \brief A mass operator
     */
    Mass myMass;


};//IP2NodeSolver

};//hfox

#endif//IP2NODESOLVER_H
