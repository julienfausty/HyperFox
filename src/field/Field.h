#ifndef FIELD_H
#define FIELD_H

#include "Mesh.h"
#include "Interpolator.h"

/*!
 * \brief interface class that allows access to the values of a field on a mesh.
 */

class Field{
  public:
    //Destructor
    virtual ~Field();
    /*!
     * \brief sets a pointer to the Mesh of the field.
     */
    virtual void setMesh(Mesh * pmesh_candidate) const=0;
    /*!
     * \brief gets a pointer to the Mesh of the field.
     */
    virtual Mesh * getMesh() const=0;
    /*!
     * \brief returns the length of the values vector.
     */
    virtual int getLength() const =0;
  protected:
    /*!
     * \brief pointer to the mesh the field is defined on.
     */
    Mesh * pmesh;
};

#endif
