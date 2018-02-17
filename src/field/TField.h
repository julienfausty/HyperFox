#ifndef TFIELD_H
#define TFIELD_H

#include <vector>
#include "Field.h"
#include "Mesh.h"
#include "Interpolator.h"

/*!
 * \brief class that allows access to the values of a field on a mesh.
 */

template <class T>
class TField : public Field{
  public:
    //Constructors
    /*!
     * \brief empty constructor
     *
     * sets the pointers to nulls.
     */
    TField();
    TField(Mesh * pmesh);
    TField(Mesh * pmesh, std::vector<T> & values_candidate);
    //Destructor
    virtual ~TField();
    /*!
     * \brief gets a pointer to the Mesh of the field.
     */
    Mesh * getMesh() const;
    /*!
     * \brief gets a const pointer to the entire field.
     */
    const std::vector<T> * getValues() const;
    /*!
     * \brief sets the entire field by copying the candidate field.
     */
    void setValues(const std::vector<T> & values_candidate);
    /*!
     * \brief get the value of the field at the ith node.
     */
    virtual T getValue(int i) const;
    /*!
     * \brief get the value of the field at the point coordinate.
     */
    virtual T interpolate(const std::vector<double> &point) const;
    /*!
     * \brief set the value of the field at the ith node.
     */
    virtual void setValue(int i, const T value_candidate);
    /*!
     * \brief returns the length of the values vector.
     */
    virtual int getLength() const;
  protected:
    /*!
     * \brief function used by the constructors.
     */
    virtual void construct(Mesh * pmesh, std::vector<T> & values_candidate);
    /*!
     * \brief storage of the field.
     */
    std::vector<T> values;
    /*!
     * \brief pointer to the mesh the field is defined on.
     */
    Mesh * pmesh;
};

#endif
