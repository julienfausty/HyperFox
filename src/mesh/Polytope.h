#ifndef POLYTOPE_H
#define POLYTOPE_H

#include <vector>
#include <iostream>
#include "ErrorHandle.h"

namespace hfox{

/*!
 * \brief an element of volume of a space delimited by facets.
 *
 * Any k-polytope as (k-1)-polytopes for facets and so on towards its vertices.
 * Basically any volume enclosed by a set of points and straight lines between 
 * these points.
 */

class Polytope{
  public:
    // Destructor
    virtual ~Polytope(){};
    // Setters
    /*!
     * \brief sets the points list.
     */
    virtual void setPoints(std::vector< std::vector<double> > & 
        points_candiate);
    /*!
     * \brief sets the faces list.
     */
    virtual void setFaces(std::vector< std::vector<int> > & faces_candidate);
    // Getters
    /*!
     * \brief gets the points list.
     */
    virtual const std::vector< std::vector<double> > * getPoints() const;
    /*!
     * \brief gets the faces list.
     */
    virtual const std::vector< std::vector<int> > * getFaces() const;
    /*!
     * \brief gets the number of points.
     */
    virtual int getNumberPoints() const;
    /*!
     * \brief gets the number of points.
     */
    virtual int getNumberFaces() const;
    /*!
     * \brief get the dimension of the polytope.
     */
    virtual int getDimension() const;
    // Calculators
    /*!
     * \brief calculate the volume enclosed by the polytope.
     * @return the volume of the polytope.
     */
    virtual double calcVolume() const =0;
    /*!
     * \brief calculate a list of the surfaces of the faces of the polytope.
     * @return a list of the surfaces of the faces of the polytope.
     */
    virtual std::vector<double> calcSurfaces() const =0;
    /*!
     * \brief calculate a list of the edge lengths of the faces of the polytope.
     * @return a list of the edge lengths of the faces of the polytope.
     */
    virtual std::vector< std::vector<double> > calcEdgeLengths() const =0;
    /*!
     * \brief calculate the total surface of the faces of the polytope.
     * @return the total surface of the faces of the polytope.
     */
    virtual double calcSurface() const =0;
    /*!
     * \brief calculate the total edge length of the polytope.
     * @return the total edge length of the polytope.
     */
    virtual double calcEdgeLength() const =0;
  protected:
    /*!
     * The list of points making the polytope.
     */
    std::vector< std::vector<double> > points;
    /*!
     * The list of faces of the polytope. Each face is the connectivity of the 
     * points by index.
     */
    std::vector < std::vector<int> > faces;
    /*!
     * dimension of the space delimited by the polytope.
     */
    int dimension;
};

} //hfox

#endif
