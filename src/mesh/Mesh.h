#ifndef MESH_H
#define MESH_H

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <moab/Core.hpp>
#include <moab/Skinner.hpp>
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
    /*! \brief A constructor for the mesh with points and connectivity.
     *
     * @param point_candidate
     * @param connectivity_candidate a pointer to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    Mesh(int dim, int order, std::string geom, moab::EntityHandle meshset_candidate);
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
    /*!
     * \brief A method for directly setting the meshset with moab compatible applications
     *
     * @param meshset_candidate the candidate for the meshset
     */
    void setMeshSet(moab::EntityHandle meshset_candidate);
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     *
     * @param points the pointer to the vector to fill with vertices
     */
    void getPoints(std::vector< std::vector<double> > * points) const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getCells(std::vector< std::vector<int> > * connectivity) const;
    /*! \brief A method for geting a const pointer to the connectivity fo the faces of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getFaces(std::vector< std::vector<int> > * faces) const;
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     *
     * @param points the pointer to the vector to fill with vertices
     */
    void getSlicePoints(const std::vector<int> & slice, std::vector< std::vector<double> > * points) const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getSliceCells(const std::vector<int> & slice, std::vector< std::vector<int> > * connectivity) const;
    /*! \brief A method for geting a const pointer to the connectivity fo the faces of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getSliceFaces(const std::vector<int> & slice, std::vector< std::vector<int> > * faces) const;
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
    /*! \brief A method for accessing an individual face.
     *
     * @param i index of the face of which one wishes to get the connectivity 
     * information. 
     * @param cell pointer to fill with node indexes
     */
    void getFace(int i, std::vector<int> * face) const;
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
    /*! \brief get the number of faces of the mesh
     */
    int getNumberFaces() const;
    /*! \brief get the dimension of the space the mesh is discretizing.
     */
    int getDimension() const;
    /*!
     * \brief get an element to face map
     */
    const std::vector< std::vector<int> > * getCell2FaceMap() const;
    /*!
     * \brief get the faces of a specific cell
     *
     * @param i index of cell
     */
    const std::vector<int> * getCell2Face(int i) const;
    /*!
     * \brief get a face to element map
     *
     * @param face2CellMap pointer to fill with cells adjoining faces
     */
    const std::vector< std::vector<int> > * getFace2CellMap() const;
    /*!
     * \brief get the cells adjoining a specific face
     *
     * @param i index of face
     */
    const std::vector<int> * getFace2Cell(int i) const;
    /*!
     * \brief get the indexes of the faces on the boundary
     */
    const std::set<int> * getBoundaryFaces() const;
  protected:
    /*!
     * \brief small method for initializing the mbInterface
     */
    void initializeMBInterface();
    /*! \brief compute and set both the inner and outer faces of the mesh.
     */
    void computeFaces();
    /*! \brief compute the face to cell map of the mesh.
     */
    void computeFace2CellMap();
    /*! \brief compute the cell to face map of the mesh.
     */
    void computeCell2FaceMap();
    /*! \brief compute the boundary of the mesh.
     */
    void computeBoundary();
    /*!
     * \brief small method for finding the internal cell index of a face
     */
    int getFaceIndex(const moab::EntityHandle & cell, const moab::EntityHandle & face) const;
    /*!
     * \brief small method for determining the moab element type
     */
    moab::EntityType determineMOABType() const;
    /*!
     * \brief small method for determining the moab face type
     */
    moab::EntityType determineMOABTypeFace() const;
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
    /*!
     * \brief the face 2 cell map
     */
    std::vector< std::vector<int> > face2CellMap;
    /*!
     * \brief the face 2 cell map
     */
    std::vector< std::vector<int> > cell2FaceMap;
    /*!
     * \brief the indexes of the boundary faces of the mesh
     */
    std::set<int> boundaryFaces;
}; //Mesh

} //hfox

#endif
