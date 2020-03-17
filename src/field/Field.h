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
    virtual int getLength() const{return nValues;};
    /*!
     * \brief allocates the memory necessary for a field of certain length.
     *
     * @param length integer with size of field
     */
    virtual int allocate(int length) =0;
  protected:
    /*!
     * \brief number of values
     */
    int nValues;
    /*!
     * \brief pointer to the mesh the field is defined on.
     */
    Mesh * pmesh;
};//Field

} //hfox

#endif//FIELD_H
