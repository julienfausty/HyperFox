#ifndef MESH_H
#define MESH_H

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <moab/Core.hpp>
#include <moab/Skinner.hpp>
#include "ReferenceElement.h"
#include "ErrorHandle.h"

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
    /*! \brief A constructor for the mesh with just parameters to the reference element.
     *
     * @param dim the dimension of the mesh
     * @param order the FE order of the mesh
     * @param geom a string with the type of element
     */
    Mesh(int dim, int order, std::string geom);
    /*! \brief A constructor for the mesh with points and connectivity.
     *
     * @param dim the dimension of the mesh
     * @param order the FE order of the mesh
     * @param geom a string with the type of element
     * @param point_candidate
     * @param connectivity_candidate a pointer to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    Mesh(int dim, int order, std::string geom, int dimPointSpace, std::vector< double > &  point_candidate,
        std::vector< int > & connectivity_candidate);
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
    void setMesh(int dimPointSpace, std::vector<double> & points_candidate, 
        std::vector<int> & connectivity_candidate);
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     */
    const std::vector<double> * getPoints() const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     */
    const std::vector<int> * getCells() const;
    /*! \brief A method for geting a const pointer to the connectivity fo the faces of the mesh.
     */
    const std::vector<int> * getFaces() const;
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     *
     * @param points the pointer to the vector to fill with vertices
     */
    void getSlicePoints(const std::vector<int> & slice, std::vector< const double * > * points) const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getSliceCells(const std::vector<int> & slice, std::vector< const int * > * connectivity) const;
    /*! \brief A method for geting a const pointer to the connectivity fo the faces of the mesh.
     *
     * @param connectivity the pointer to the vector to fill with connectivity
     */
    void getSliceFaces(const std::vector<int> & slice, std::vector< const int * > * faces) const;
    /*! \brief A method for accessing a specific node of the mesh.
     *
     * @param i index of the point one wishes to acquire.
     */
    const double * getPoint(int i) const;
    /*! \brief A method for accessing an individual cell.
     *
     * @param i index of the cell of which one wishes to get the connectivity 
     * information. 
     * @param cell pointer to fill with node indexes
     */
    const int * getCell(int i) const;
    /*! \brief A method for accessing an individual face.
     *
     * @param i index of the face of which one wishes to get the connectivity 
     * information. 
     */
    const int * getFace(int i) const;
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
    /*! \brief get the dimension of the space the mesh is embedded in.
     */
    int getNodeSpaceDimension() const;
    /*! \brief get the number of faces per cell.
     */
    int getNumFacesPerCell() const;
    /*!
     * \brief get an element to face map
     */
    const std::vector<int> * getCell2FaceMap() const;
    /*!
     * \brief get the faces of a specific cell
     *
     * @param i index of cell
     */
    const int * getCell2Face(int i) const;
    /*!
     * \brief get a face to element map
     */
    const std::vector<int> * getFace2CellMap() const;
    /*!
     * \brief get the cells adjoining a specific face
     *
     * @param i index of face
     */
    const int * getFace2Cell(int i) const;
    /*!
     * \brief get the indexes of the faces on the boundary
     */
    const std::set<int> * getBoundaryFaces() const;
  protected:
    /*! \brief compute and set both the inner and outer faces of the mesh.
     */
    void computeFaces();
    /*! \brief compute the face to cell map of the mesh.
     */
    void computeFace2CellMap(moab::Interface * mbInterface, moab::EntityHandle & meshset);
    /*! \brief compute the cell to face map of the mesh.
     */
    void computeCell2FaceMap(moab::Interface * mbInterface, moab::EntityHandle & meshset, ReferenceElement & skelRefEl);
    /*! \brief compute the boundary of the mesh.
     */
    void computeBoundary(moab::Interface * mbInterface, moab::EntityHandle & meshset);
    /*!
     * \brief small method for finding the internal cell index of a face
     */
    int getFaceIndex(moab::Interface * mbInterface, 
        const moab::EntityHandle & cell, const moab::EntityHandle & face, ReferenceElement & skelRefEl) const;
    /*!
     * \brief small method for determining the moab element type
     */
    moab::EntityType determineMOABType() const;
    /*!
     * \brief small method for determining the moab face type
     */
    moab::EntityType determineMOABTypeFace() const;
    /*!
     * \brief where the node coordinates are stored
     */
    std::vector<double> nodes;
    /*!
     * \brief dimension of the space of the nodes
     */
    int dimNodeSpace;
    /*!
     * \brief the number of nodes
     */
    int nNodes;
    /*!
     * \brief where the connectivity of the cells are stored
     */
    std::vector<int> cells;
    /*!
     * \brief the number of nodes per cell
     */
    int nNodesPerCell;
    /*!
     * \brief the number of cells
     */
    int nCells;
    /*!
     * \brief where the connectivity of the faces are stored
     */
    std::vector<int> faces;
    /*!
     * \brief the number of nodes per face
     */
    int nNodesPerFace;
    /*!
     * \brief the number of faces
     */
    int nFaces;
    /*!
     * \brief the reference element of the mesh
     */
    ReferenceElement * refElement;
    /*!
     * \brief the face 2 cell map
     */
    std::vector<int> face2CellMap;
    /*!
     * \brief the cell 2 face map
     */
    std::vector<int> cell2FaceMap;
    /*!
     * \brief number of faces per cell
     */
    int nFacesPerCell;
    /*!
     * \brief the indexes of the boundary faces of the mesh
     */
    std::set<int> boundaryFaces;
}; //Mesh

} //hfox

#endif
