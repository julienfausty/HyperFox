#ifndef MODEL_H
#define MODEL_H

#include <map>
#include <string>
#include "ReferenceElement.h"
#include "DenseEigen.h"


namespace hfox{

/*!
 * \brief An interface class that describes the necessary elements of a model
 *
 * A Model is supposed to be an object that represents the system one is trying to simulate. As such, it must 
 * implement some element level mechanics of the simulation. It should be useful regardless of the FE method or
 * assembly.
 */

class Model{
  public:
    /*!
     * \brief a method to set the local fields
     *
     * @param pointer to a map of names and corresponding local values of fields
     */
    virtual void setFieldMap(const std::map<std::string, std::vector<double> > * fm){fieldMap = fm;};
    /*!
     * \brief a method to set the reference element
     */
    virtual void setReferenceElement(const ReferenceElement * re){refEl = re;};
    /*!
     * \brief a method to set the current element nodes
     */
    virtual void setElementNodes(const std::vector< std::vector<double> > * ens){elementNodes = ens;};
    /*!
     * \brief a method where all the relevant values must be computed
     */
    virtual void compute() =0;
  protected:
    /*!
     * \brief the current element nodes
     */
    const std::vector< std::vector<double> > * elementNodes;
    /*!
     * \brief the reference element of the FE mesh
     */
    const ReferenceElement * refEl;
    /*!
     * \brief a pointer to the map containing all the local fields
     */
    const std::map<std::string, std::vector<double> > * fieldMap;

}; //Model

}//hfox

#endif //MODEL_H
