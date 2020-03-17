#ifndef REFERENCEELEMENT_H
#define REFERENCEELEMENT_H

namespace hfox{

/*!
 * \brief The main reference element class
 *
 * The basis of the finite element method. The class can represent higher order elements
 * with possibly different geometries in arbitrary dimensions.
 */

class ReferenceElement{

  public:
  protected:
    /*!
     * \brief the space dimension of the element
     */
    int dimension;
    /*!
     * \brief the polynomial order of the element
     */
    int order;
    /*!
     * \brief the type of element
     */
    int order;
}//ReferenceElement

}//hfox

#endif//REFERENCEELEMENT_H
