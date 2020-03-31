#ifndef MESH_H
#define MESH_H

#include <vector>
#include <array>
#include <iostream>
#include <moab/Core.hpp>
#include "ReferenceElement.h"

namespace hfox{

/*!
 * \brief The Mesh class.
 *
 * Any version of variable discretization must be a mesh object (space or 
 * time). This class is represents that discretization. 
 */

class Mesh{
  public:
    //Constructors
    /*! \brief An empty constructor for the mesh.
     */
    Mesh();
    /*! \brief A constructor for the mesh with points and connectivity.
     *
     * @param point_candidate
     * @param connectivity_candidate a pointer to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    Mesh(int dim, int order, std::string geom, std::vector< std::vector< double > > &  point_candidate,
        std::vector< std::vector< int > > & connectivity_candidate);
    //Destructors
    /*! \brief The destructor for the mesh.
     */
    ~Mesh();
    /*!
     * \brief A method for setting the reference element.
     *
     * @param dim the dimension of the space
     * @param order the order of the polynomial interpolation
     * @param geom the geometry of the elements
     */
    void setReferenceElement(int dim, int order, std::string geom);
    /*! \brief A method for setting the points and connectivity of the mesh.
     *
     * @param point_candidate a reference to a vector of coordinates of dimension
     * d.
     * @param connectivity_candidate a reference to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    void setMesh(std::vector< std::vector<double> > & points_candidate, 
        std::vector< std::vector<int> > & connectivity_candidate);
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     *
     * @param points the pointer to the vector to fill with vertices
     */
    void getPoints(std::vector< std::vector<double> > * points) const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getConnectivity(std::vector< std::vector<int> > * connectivity) const;
    /*! \brief A method for accessing a specific node of the mesh.
     *
     * @param i index of the point one wishes to acquire.
     * @param point pointer to fill with coordinates
     */
    void getPoint(int i, std::vector<double> * point) const;
    /*! \brief A method for accessing an individual cell.
     *
     * @param i index of the cell of which one wishes to get the connectivity 
     * information. 
     * @param cell pointer to fill with node indexes
     */
    void getCell(int i, std::vector<int> * cell) const;
    /*!
     * \brief get a const pointer to the reference element
     */
    const ReferenceElement * getReferenceElement() const;
    // Helper functions
    /*! \brief get the number of nodes of the mesh
     */
    int getNumberPoints() const;
    /*! \brief get the number of cells of the mesh
     */
    int getNumberCells() const;
    /*! \brief get the dimension of the space the mesh is discretizing.
     */
    int getDimension() const;
    /*! \brief A method for geting a const pointer to the inner faces of the mesh.
     *
     * @param innerFaces pointer to a vector to fill with inner faces
     */
    void getInnerFaces(std::vector< std::array<int, 4> > * innerFaces) const;
    /*! \brief A method for geting a const pointer to an inner face of the mesh.
     *
     * @param i index of face
     * @param innerFace pointer to fill with inner face
     */
    void getInnerFace(int i, std::array<int, 4> * innerFace) const;
    /*! \brief A method for geting a const pointer to the outer faces of the mesh.
     *
     * @param outerFaces pointer to a vector to fill with outer faces
     */
    void getOuterFaces(std::vector< std::array<int, 2> > * outerFaces) const;
    /*! \brief A method for geting a const pointer to an outer face of the mesh.
     *
     *  @param i index to outer face
     *  @param outerFace a pointer to fill with outer face info
     */
    void getOuterFace(int i, std::array<int, 2> * outerFace) const;
  protected:
    /*! \brief compute and set both the inner and outer faces of the mesh.
     */
    void computeFaces();
    /*!
     * \brief small method for initializing the mbInterface
     */
    void initializeMBInterface();
    /*!
     * \brief small method for determining the moab element type
     */
    moab::EntityType determineMOABType() const;
    /*!
     * \brief the interface to the moab description
     */
    moab::Interface * mbInterface;
    /*!
     * \brief moab entity handle to the mesh set
     */
    moab::EntityHandle meshset;
    /*!
     * \brief the reference element of the mesh
     */
    ReferenceElement * refElement;
}; //Mesh

} //hfox

#endif
