#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <vector>

namespace hfox{

/*!
 * \brief a helper class for computing convex hulls in arbitrary dimensions.
 *
 * Should be utilized as a toolbox for computing a convex hull given a set of 
 * points.
 */

template<int d> class Simplex;

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
     * \brief for setting vertixes pointer.
     */
    void setFaces(std::vector< std::vector<int> > * pfaces_cand);
    /*!
     * \brief for setting the algorithm type to use.
     */
    void setAlgoType(CHullAlgo typeAlgo_cand);
    //Getters
    /*!
     * \brief for getting the vertixes pointer.
     */
    const std::vector< std::vector<double> > * getVertexes() const;
    /*!
     * \brief for getting vertixes pointer.
     */
    std::vector< std::vector<int> > * getFaces() const;
    /*!
     * \brief for getting the algorithm type to use.
     */
    CHullAlgo getAlgoType() const;
    // Calculators
    /*!
     * \brief the entry function to the convex hull computation.
     */
    void computeConvexHull();
    /*!
     * \brief an implementation of the QuickHull algorithm (Barber et al., 1996)
     */
    void quickHull();
  protected:
    int dimension;
    const std::vector< std::vector<double> > * pvertexes;
    std::vector< std::vector<int> > * pfaces;
    CHullAlgo typeAlgo;
};

} //hfox

#endif
