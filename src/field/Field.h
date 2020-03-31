#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include "Mesh.h"

namespace hfox{

/*!
 * \brief Interface class to a field on a mesh.
 */

class Field {
  public:
    //Destructor
    virtual ~Field();
    /*!
     * \brief gets a pointer to the Mesh of the field.
     */
    Mesh * getMesh() const{return pmesh;};
    /*!
     * \brief gets a pointer to the Mesh of the field.
     */
    void setMesh(Mesh * m){pmesh = m;};
    /*!
     * \brief returns the length of the values vector.
     */
    int getLength() const{return values.size();};
    /*!
     * \brief allocates the memory necessary for a field of certain length.
     *
     * @param length integer with size of field
     */
    int allocate(int length){values.resize(length);};
    /*!
     * \brief get pointer to the value list
     */
    std::vector<double> * getValues(){return &values;};
  protected:
    /*!
     * \brief list of values of the field
     */
    std::vector<double> values;
    /*!
     * \brief pointer to the mesh the field is defined on.
     */
    Mesh * pmesh;
};//Field

} //hfox

#endif//FIELD_H
