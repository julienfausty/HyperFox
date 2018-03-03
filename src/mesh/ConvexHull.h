#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <vector>
#include "ErrorHandle.h"
#include "Simplex.h"
#include "DenseEigen.h"

namespace hfox{

/*!
 * \brief a helper class for computing convex hulls in arbitrary dimensions.
 *
 * Should be utilized as a toolbox for computing a convex hull given a set of 
 * points.
 */

enum CHullAlgo{ QuickHull };

class ConvexHull{
  public:
    // Constructor
    ConvexHull() : typeAlgo(QuickHull){};
    //Setters
    /*!
     * \brief for setting vertixes pointer.
     */
    void setVertexes(const std::vector< std::vector<double> > * pvertixes_cand);
    /*!
     * \brief for setting faces pointer.
     */
    void setFaces(std::vector< std::vector<int> > * pfaces_cand);
    /*!
     * \brief for setting the algorithm type to use.
     */
    void setAlgoType(CHullAlgo typeAlgo_cand);
    //Getters
    /*!
     * \brief for getting the algorithm type to use.
     */
    CHullAlgo getAlgoType() const;
    /*!
     * \brief for getting the algorithm type to use.
     */
    int getDimension() const;
    // Calculators
    /*!
     * \brief the entry function to the convex hull computation.
     */
    void computeConvexHull();
    /*!
     * \brief calculate the signed distance from face.
     * @param face indexes of the points comprising the face.
     * @param point the point under consideration.
     * @return the signed distance from the point to the face with regards to its normal.
     */
    double calcDistFace(std::vector<int> & face, 
        std::vector<double> & point) const;
    /*!
     * \brief calculate the oriented normal to the face.
     * @param face indexes of the points comprising the face.
     * @return the oriented unit normal to the face.
     */
    EVector calcNormalVector(std::vector<int> & face) const;
    /*!
     * \brief remove the column i from the matrix mat and return a new matrix.
     */
    EMatrix removeColumn(const EMatrix & mat, unsigned int colToRemove) const;
  protected:
    /*!
     * \brief an implementation of the QuickHull algorithm (Barber et al., 1996)
     */
    void quickHull();
    /*!
     * dimension of the space.
     */
    int dimension;
    /*!
     * pointer to the vertices.
     */
    const std::vector< std::vector<double> > * pvertexes;
    /*!
     * pointer to the faces supposed to be computed by the class.
     */
    std::vector< std::vector<int> > * pfaces;
    /*!
     * enumerator type for supported algorithms.
     */
    CHullAlgo typeAlgo;
};

} //hfox

#endif
