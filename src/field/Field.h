#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include "Mesh.h"
#include "FieldTypes.h"

namespace hfox{

/*!
 * \brief Interface class to a field on a mesh.
 */

class Field {
  public:
    /*!
     * \brief empty constructor
     */
    Field();
    /*!
     * \brief mesh constructor
     *
     * @param mesh pointer to the mesh
     */
    Field(Mesh * mesh);
    /*!
     * \brief allocating conatructor
     *
     * @param mesh pointer to the mesh
     * @param t the type of field (defined on nodes, edges, faces or cells)
     * @param nObjPerEnt number of objects per geometry entity
     * @param nValPerObj number of values per each object
     */
    Field(Mesh * mesh, FieldType t, int nObjPerEnt, int nValPerObj);
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
     */
    void allocate();
    /*!
     * \brief get pointer to the value list
     */
    std::vector<double> * getValues(){return &values;};
    /*!
     * \brief get pointer to the number of objects per entity
     */
    FieldType * getFieldType(){return &type;};
    /*!
     * \brief get pointer to the number of entities
     */
    int * getNumEntities(){return &numEntities;};
    /*!
     * \brief get pointer to the number of objects per entity
     */
    int * getNumObjPerEnt(){return &numObjPerEnt;};
    /*!
     * \brief get pointer to the number of values per object
     */
    int * getNumValsPerObj(){return &numValsPerObj;};
    /*!
     * \brief compute the number of entities from the field type
     */
    void computeNumEntities();
  protected:
    /*!
     * \brief list of values of the field
     */
    std::vector<double> values;
    /*!
     * \brief pointer to the mesh the field is defined on.
     */
    Mesh * pmesh;
    /*!
     * \brief the type of field
     */
    FieldType type;
    /*!
     * \brief the number of objects per entity
     */
    int numEntities;
    /*!
     * \brief the number of objects per entity
     */
    int numObjPerEnt;
    /*!
     * \brief the number of values per object
     */
    int numValsPerObj;
};//Field

} //hfox

#endif//FIELD_H
