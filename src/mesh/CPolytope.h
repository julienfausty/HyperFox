#ifndef CPOLYTOPE_H
#define CPOLYTOPE_H

#include <vector>
#include <iostream>
#include "Polytope.h"
#include "ErrorHandle.h"

namespace hfox{

// Forward declaration.
class ConvexHull;

/*!
 * \brief a convex element of volume of a space delimited by facets.
 *
 * Any k-polytope as (k-1)-polytopes for facets and so on towards its vertices.
 * Basically any CONVEX volume enclosed by a set of points and straight lines 
 * between these points.
 */

class CPolytope : public Polytope{
  public:
    CPolytope();
    CPolytope(std::vector< std::vector<double> > & points_cand);
    // Setters
    virtual void setPoints(std::vector< std::vector<double> > & 
        points_candidate){Polytope::setPoints(points_candidate);};
    /*!
     * \brief should never be called for CPolytopes but calculated.
     */
    virtual void setFaces(std::vector< std::vector<int> > & faces_candidate);
    /*!
     * \brief sets the faces list by calculating the ConvexHull.
     */
    virtual void setFaces();
    // Calculators
    /*!
     * \brief calculate the convex hull of the polytope in faces. 
     */
    virtual void calcConvexHull();
    /*!
     * \brief calculate the signed distance from face (- inside, + outside)
     * @param index_face the index of the considered face in the faces list.
     * @param point the point under consideration.
     * @return the distance from the point to the face.
     */
    virtual double calcDistFace(int index_face, 
        std::vector<double> & point) const;
    /*!
     * \brief calculate the volume enclosed by the polytope.
     * @return the volume of the polytope.
     */
    virtual double calcVolume() const;
    /*!
     * \brief calculate a list of the surfaces of the faces of the polytope.
     * @return a list of the surfaces of the faces of the polytope.
     */
    virtual std::vector<double> calcSurfaces() const;
    /*!
     * \brief calculate a list of the edge lengths of the faces of the polytope.
     * @return a list of the edge lengths of the faces of the polytope.
     */
    virtual std::vector< std::vector<double> > calcEdgeLengths() const;
    /*!
     * \brief calculate the total surface of the faces of the polytope.
     * @return the total surface of the faces of the polytope.
     */
    virtual double calcSurface() const;
    /*!
     * \brief calculate the total edge length of the polytope.
     * @return the total edge length of the polytope.
     */
    virtual double calcEdgeLength() const;
    /*!
     */
    //Checkers
    /*!
     * \brief check if point is within the bounds of the polytope.
     */
    virtual bool isPointInside(const std::vector<double> & point) const;
};

} //hfox

#endif //CPOLYTOPE_H
