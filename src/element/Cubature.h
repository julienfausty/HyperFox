#ifndef CUBATURE_H
#define CUBATURE_H

#include <string>
#include <vector>
#include "ElementGeometry.h"

namespace hfox{

/*!
 * \brief Class for curbature rules for polynomials
 *
 * A class that computes curbature rules for arbitrary degree polynomials in arbitrary dimensions.
 * These rules are taken from the computations performed in 
 *
 * F.D. Witherden , P.E. Vincent, "On the identification of symmetric quadrature rules for finite 
 * element methods", Computers and Mathematics with Applications 69 (2015) 1232–1241
 */


class Cubature{
  public:
    //Constructors
    /*!
     * \brief the default constructor for the reference element
     */
    Cubature(int dim, int ord, elementGeometry geom);
    //Destructor
    ~Cubature();
    //Exposed methods
    /*!
     * \brief get dimension of space
     */
    int getDimension() const;
    /*!
     * \brief get polynomial interpolation order
     */
    int getOrder() const;
    /*!
     * \brief get element geometry
     */
    elementGeometry getGeometry() const;
    /*!
     * \brief get number of IPs
     */
    int getNumIPs() const;
    /*!
     * \brief get a const pointer to the IP coordinates
     */
    const std::vector< std::vector<double> > * getIPCoords() const;
    /*!
     * \brief get a const pointer to the IP weights
     */
    const std::vector<double> * getIPWeights() const;
  protected:
    /*!
     * \brief method for determining the number of integration points
     */
    void determineNumIPs();
    /*!
     * \brief method for determining the coordinates of the IPs
     */
    void determineIPs();
    /*!
     * \brief method for determining the weights of the IPs
     */
    void determineWeightsIPs();
    /*!
     * \brief the space dimension of the element
     */
    int dimension;
    /*!
     * \brief the polynomial order of the element
     */
    int order;
    /*!
     * \brief the geometry of the element
     */
    elementGeometry geometry;
    /*!
     * \brief the number of integration points
     */
    int nIPs;
    /*!
     * \brief coordinates of the integration points
     */
    std::vector< std::vector<double> > ipCoords;
    /*!
     * \brief weights of the integration points
     */
    std::vector<double> ipWeights;
}; //Cubature

}//hfox

#endif //CUBATURE_H
