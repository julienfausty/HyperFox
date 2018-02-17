#ifndef OPERATOR_H
#define OPERATOR_H

#include "Mesh.h"

/*!
 * \brief an object for performing calculations on fields.
 *
 * The operator operates on a field or multiple fields to compute a value or 
 * another field.
 */

class Operator{
  public:
    // Destructor
    virtual ~Operator();
  protected:
    const Mesh * pmesh;
}

#endif
