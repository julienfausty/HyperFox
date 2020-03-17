#ifndef MESH_H
#define MESH_H

#include <vector>
#include <array>
#include "Io.h"
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
    Mesh(std::vector< std::vector< double > > &  point_candidate,
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
    void setReferenceElement(int dim, int order, elementGeometry geom);
    /*! \brief A method for setting the points of the mesh.
     *
     * @param point_candidate a reference to a vector of coordinates of dimension
     * d.
     */
    void setPoints(std::vector< std::vector<double> > &
        points_candidate);
    /*! \brief A method for setting the connectivity of the mesh.
     *
     * @param connectivity_candidate a reference to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    void setConnectivity(std::vector< std::vector<int> > &
        connectivity_candidate);
    /*! \brief A method for geting a const pointer to the nodes of the mesh.
     */
    const std::vector< std::vector<double> > * getPoints() const;
    /*! \brief A method for geting a const pointer to the connectivity of the mesh.
     */
    const std::vector< std::vector<int> > * getConnectivity() const;
    /*! \brief A method for accessing a specific node of the mesh.
     *
     * @param i index of the point one wishes to acquire. 
     */
    const std::vector<double> * getPoint(int i) const;
    /*! \brief A method for accessing an individual cell.
     *
     * @param i index of the cell of which one wishes to get the connectivity 
     * information. 
     */
    const std::vector<int> * getCell(int i) const;
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
    /*! \brief A method for geting a const pointer to the faces of the mesh.
     *
     *  Warning: if the faces of the mesh have not been computed, this method
     *  will compute them!
     */
    const std::vector< std::array<int, 4> > * getFaces();
    /*! \brief A method for geting a const pointer to a face of the mesh.
     *
     *  Warning: if the faces of the mesh have not been computed, this method
     *  will compute them!
     */
    const std::array<int, 4> * getFace(int i);
  protected:
    /*! \brief get the dimension of the space the mesh is discretizing.
     */
    void computeFaces();
    /*!
     * The vector of coordinates describing the nodes of the mesh.
     */
    std::vector< std::vector<double> > points;
    /*!
     * The vector of vectors of indexes of points holding the 
     * connectivity information of the mesh.
     */
    std::vector< std::vector<int> > connectivity;
    /*!
     * \brief the reference element of the mesh
     */
    ReferenceElement refElement;
    /*!
     * The vector of integers identifying the points on the boundary
     */
    std::vector<int> boundaryPoints;
    /*!
     * The vector of arrays of indexes of cells describing the inner faces. The structure is
     * [ el1 faceofel1 el2 faceofel2 ]
     */
    std::vector< std::array<int, 4> > innerFaces;
    /*!
     * The vector of arrays of indexes of cells describing the outer faces. The structure is
     * [ el1 faceofel1 ]
     */
    std::vector< std::array<int, 2> > outerFaces;
    /*!
     * The dimension of the space.
     */
    int dimension;
    /*!
     * A boolean indicating if the faces of the mesh have been computed.
     */
    bool facesExist;
}; //Mesh

} //hfox

#endif
