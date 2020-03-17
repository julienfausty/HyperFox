#ifndef MODEL_H
#define MODEL_H

#include <map>
#include <string>
#include "ReferenceElement.h"
#include "Field.h"
#include "DenseEigen.h"


namespace hfox{

/*!
 * \brief An interface class that describes the necessary elements of a model
 *
 * A Model is supposed to be an object that represents the system one is trying to simulate. As such, it must 
 * implement the element level mechanics of the simulation. It should be useful regardless of the FE method or
 * assembly.
 */

class Model{
  public:
    /*!
     * \brief a method to set the fields
     */
    virtual void setFieldMap(std::map<std::string, Field * > * fm){fieldMap = fm;};
    /*!
     * \brief a method to set the reference element
     */
    virtual void setReferenceElement(const ReferenceElement * re){refEl = re;};
    /*!
     * \brief a method to set the current element nodes
     */
    virtual void setElementNodes(std::vector< std::vector<double> * > & ens){elementNodes = &ens;};
    /*!
     * \brief compute the local matrix
     */
    virtual void computeLocalMatrix()=0;
    /*!
     * \brief compute the local right hand side
     */
    virtual void computeLocalRHS()=0;
    /*!
     * \brief compute the local matrix
     */
    virtual const EMatrix getLocalMatrix() const{return localMatrix;};
    /*!
     * \brief compute the local right hand side
     */
    virtual const EVector getLocalRHS() const{return localRHS;};
  protected:
    /*!
     * \brief the current element nodes
     */
    std::vector< std::vector<double> * > * elementNodes;
    /*!
     * \brief the reference element of the FE mesh
     */
    const ReferenceElement * refEl;
    /*!
     * \brief a pointer to the map containing all the fields
     */
    std::map<std::string, Field * > * fieldMap;
    /*!
     * \brief the local element matrix
     */
    EMatrix localMatrix;
    /*!
     * \brief the local element right hand side
     */
    EVector localRHS;


}; //Model

}//hfox

#endif //MODEL_H
