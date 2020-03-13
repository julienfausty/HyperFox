#ifndef MESH_H
#define MESH_H

#include <vector>
#include "Io.h"

namespace hfox{

/*!
 * \brief The Mesh class.
 *
 * Any version of variable discretization must be a mesh object (space or 
 * time). This class is represents that discretization. The d 
 * variable is the dimension of the discretization.
 */

class Mesh{
  public:
    //Constructors
    /*! \brief An empty constructor for the mesh.
     */
    Mesh(){ ; };
    /*! \brief A constructor for the mesh that takes an Input/Output object as an entry.
     *
     * @param io an input/output object.
     */
    Mesh(Io io);
    /*! \brief A constructor for the mesh with points and connectivity.
     *
     * @param point_candidate
     * @param connectivity_candidate a pointer to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    Mesh(std::vector< std::vector< double > > &  ppoint_candidate,
        std::vector< std::vector< int > > & connectivity_candidate);
    //Destructors
    /*! \brief The destructor for the mesh.
     */
    ~Mesh();
    /*! \brief A method for setting the points of the mesh.
     *
     * @param point_candidate a reference to a vector of coordinates of dimension
     * d.
     */
    virtual void setPoints(std::vector< std::vector<double> > &
        points_candidate);
    /*! \brief A method for setting the connectivity of the mesh.
     *
     * @param connectivity_candidate a reference to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    virtual void setConnectivity(std::vector< std::vector<int> > &
        connectivity_candidate);
    /*! \brief A method for geting all the nodes of the mesh.
     */
    virtual const std::vector< std::vector<double> > * getPoints() const;
    /*! \brief A method for geting the connectivity of the mesh.
     */
    virtual const std::vector< std::vector<int> > * getConnectivity() const;
    /*! \brief A method for accessing a specific node of the mesh.
     *
     * @param i index of the point one wishes to acquire. 
     */
    virtual const std::vector<double> * getPoint(int i) const;
    /*! \brief A method for accessing the connectivity of a node.
     *
     * @param i index of the point of which one wishes to get the connectivity 
     * information. 
     */
    virtual const std::vector<int> * getPointConnectivity(int i) const;
    // Helper functions
    virtual int getNumberPoints() const;
    virtual int getNumberUnits() const;
    virtual int getDimension() const;
  protected:
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
     * The dimension of the space.
     */
    int dimension;
}; //Mesh

} //hfox

#endif
