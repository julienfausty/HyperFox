#ifndef REFERENCEELEMENT_H
#define REFERENCEELEMENT_H

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <numeric>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include "ElementGeometry.h"
#include "Cubature.h"
#include "ErrorHandle.h"
#include "DenseEigen.h"

namespace hfox{

/*!
 * \brief The main reference element class
 *
 * The basis of the finite element method. The class can represent higher order elements
 * with possibly different geometries in arbitrary dimensions.
 */

class ReferenceElement{

  public:
    //Constructors
    /*!
     * \brief the default constructor for the reference element
     */
    ReferenceElement(int dim, int ord, std::string geom);
    //Destructor
    ~ReferenceElement();
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
     * \brief get number of nodes
     */
    int getNumNodes() const;
    /*!
     * \brief get number of IPs
     */
    int getNumIPs() const;
    /*!
     * \brief get number of faces
     */
    int getNumFaces() const;
    /*!
     * \brief get a const pointer to the nodal coordinates
     */
    const std::vector< std::vector<double> > * getNodes() const;
    /*!
     * \brief get a const pointer to the positions of the face nodes by face in the nodes list
     */
    const std::vector< std::vector<int> > * getFaceNodes() const;
    /*!
     * \brief get a const pointer to the positions of the inner nodes in the nodes list
     */
    const std::vector<int> * getInnerNodes() const;
    /*!
     * \brief get a const pointer to the IP coordinates
     */
    const std::vector< std::vector<double> > * getIPCoords() const;
    /*!
     * \brief get a const pointer to the IP weights
     */
    const std::vector<double> * getIPWeights() const;
    /*!
     * \brief get a const pointer to the IP shape function values (inner vector are the shape funcs and outer vector is IPs)
     */
    const std::vector< std::vector<double> > * getIPShapeFunctions() const;
    /*!
     * \brief get a const pointer to the IP derivative shape function values (inner vectors are the shape funcs and outer vector is IPs)
     */
    const std::vector< std::vector< std::vector<double> > > * getIPDerivShapeFunctions() const;
    /*!
     * \brief get a const pointer to the nodal derivative shape function values (inner vectors are the shape funcs and outer vector are nodes)
     */
    const std::vector< std::vector< std::vector<double> > > * getDerivShapeFunctions() const;
    /*!
     * \brief get a const pointer to the face element
     */
    const ReferenceElement * getFaceElement() const;
    /*!
     * \brief a method for interpolating shape functions at a given point
     *
     * @param point reference to a vector of coordinates in reference space.
     */
    std::vector<double> interpolate(const std::vector<double> & point) const;
    /*!
     * \brief a method for interpolating derivatives of shape functions at a given point
     *
     * @param point reference to a vector of coordinates in reference space.
     * @param degree degree of derivative
     */
    std::vector< std::vector<double> > interpolateDeriv(const std::vector<double> & point, int degree) const;
  protected:
    /*!
     * \brief method for checking and setting element geometry
     */
    void setGeometry(std::string geom);
    /*!
     * \brief method for checking and setting dimension
     */
    void setDim(int dim);
    /*!
     * \brief method for checking and setting polynomial interpolaton order
     */
    void setOrder(int order);
    /*!
     * \brief method for determining the number of nodes
     */
    void determineNumNodes();
    /*!
     * \brief method for determining the coordinates of the nodes
     */
    void determineNodes();
    /*!
     * \brief method for determining the cubature rule in the element
     */
    void determineCubature();
    /*!
     * \brief method for determining the number of faces
     */
    void determineNumFaces();
    /*!
     * \brief method for determining the face nodes (and by extension the inner nodes)
     */
    void determineFaceNodes();
    /*!
     * \brief method for initializing face element
     */
    void setFaceElement();
    /*!
     * \brief method for computing the inverse Vandermonde matrix
     */
    void computeInverseVandermonde();
    /*!
     * \brief method for computing the shape functions at the IPs
     */
    void computeIPShapeFunctions();
    /*!
     * \brief method for computing the derivative of shape functions at the IPs
     */
    void computeIPDerivShapeFunctions();
    /*!
     * \brief method for computing the derivative of shape functions at the nodal coordinates
     */
    void computeDerivShapeFunctions();
    /*!
     * \brief a helper function for recursively generating all combinations
     */
    std::vector< std::vector<int> > generateCombinationsWithRepetition(std::vector<int> & set, int size) const;
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
     * \brief the number of nodes
     */
    int nNodes;
    /*!
     * \brief the number of faces
     */
    int nFaces;
    /*!
     * \brief a pointer to the curbature object for the element
     */
    Cubature * cubatureRule;
    /*!
     * \brief the nodal coordinates
     */
    std::vector< std::vector<double> > nodes;
    /*!
     * \brief the face nodes indexes by face
     */
    std::vector< std::vector<int> > faceNodes;
    /*!
     * \brief the indexes of the inner nodes
     */
    std::vector<int> innerNodes;
    /*!
     * \brief values of shape functions at IPs (outer vector is shape funcs and inner is IPs).
     */
    std::vector< std::vector<double> > ipShapeFunctions;
    /*!
     * \brief derivatives of shape functions at IPs (outer vector is shape funcs, middle is 
     * IPs and inner is dimensions).
     */
    std::vector< std::vector< std::vector<double> > > ipDerivShapeFunctions;
    /*!
     * \brief derivatives of shape functions at nodes (outer vector is shape funcs, middle is 
     * nodes and inner is dimensions).
     */
    std::vector< std::vector< std::vector<double> > > derivShapeFunctions;
    /*!
     * \brief reference element of a face 
     */
    ReferenceElement * faceElement;
    /*!
     * \brief matrix for computing modal basis from nodal basis
     *
     * invV_i^j = inv(P_i(x^j)) (inverse Vandermonde matrix)
     */
    EMatrix invV;
    /*!
     * \brief a function that initializes the database
     */
    static void initializeDatabase();
    /*!
     * \brief a boolean tracking the completeness of the database
     */
    static bool databaseExists;
    /*!
     * \brief the maximum allowed spatial dimension
     */
    static int maxDim;
    /*!
     * \brief the maximum allowed orders (maxOrderMap[(dim, geometry)] = maxOrder)
     */
    static std::map< std::tuple<int, elementGeometry>, int> maxOrderMap;
    /*!
     * \brief the database of element nodes (nodeDatabase[(dim, numNodes, geometry)] = nodes)
     */
    static std::map< std::tuple<int, int, elementGeometry>, std::vector< std::vector<double> > > nodeDatabase;

};//ReferenceElement

}//hfox

#endif//REFERENCEELEMENT_H
