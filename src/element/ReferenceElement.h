#ifndef REFERENCEELEMENT_H
#define REFERENCEELEMENT_H

#include <vector>
#include <string>

namespace hfox{

enum elementGeometry{
  simplex,
  orthotope
};

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
     * \brief get a const pointer to the IP coordinates
     */
    const std::vector< std::vector<double> > * getIPCoords() const;
    /*!
     * \brief get a const pointer to the IP weights
     */
    const std::vector<double> * getIPWeights() const;
    /*!
     * \brief get a const pointer to the IP shape function values
     */
    const std::vector< std::vector<double> > * getIPShapeFunctions() const;
    /*!
     * \brief get a const pointer to the IP derivative shape function values
     */
    const std::vector< std::vector< std::vector<double> > > * getIPDerivShapeFunctions() const;
    /*!
     * \brief get a const pointer to the nodal derivative shape function values
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
    std::vector<double> interpolate(std::vector<double> & point) const;
    /*!
     * \brief a method for interpolating derivatives of shape functions at a given point
     *
     * @param point reference to a vector of coordinates in reference space.
     * @param degree degree of derivative
     */
    std::vector<double> interpolateDeriv(std::vector<double> & point, int degree) const;
  protected:
    /*!
     * \brief method for determining the number of nodes
     */
    void determineNumNodes();
    /*!
     * \brief method for determining the number of integration points
     */
    void determineNumIPs();
    /*!
     * \brief method for determining the coordinates of the nodes
     */
    void determineNodes();
    /*!
     * \brief method for determining the coordinates of the IPs
     */
    void determineIPs();
    /*!
     * \brief method for determining the weights of the IPs
     */
    void determineWeightsIPs();
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
     * \brief the number of integration points
     */
    int nIPs;
    /*!
     * \brief the number of faces
     */
    int nFaces;
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
     * \brief coordinates of the integration points
     */
    std::vector< std::vector<double> > ipCoords;
    /*!
     * \brief weights of the integration points
     */
    std::vector<double> ipWeights;
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

};//ReferenceElement

}//hfox

#endif//REFERENCEELEMENT_H
