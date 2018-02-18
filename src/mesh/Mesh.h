#ifndef MESH_H
#define MESH_H

#include <vector>
#include "Polytope.h"

namespace hfox{

/*!
 * \brief The interface class for a mesh.
 *
 * Any version of variable discretization must be a mesh object (space or 
 * time). This class is the interface that each mesh object must follow. The d 
 * variable is the dimension of the discretization.
 */

class Mesh{
  public:
    /*! \brief A method for setting the points of the mesh.
     *
     * @param ppoint_candidate a pointer to a vector of coordinates of dimension
     * d.
     */
    virtual void setPoints(std::vector< std::vector<double> > &
        points_candidate)=0;
    /*! \brief A method for setting the connectivity of the mesh.
     *
     * @param ppoint_candidate a pointer to a vector of vectors holding the 
     * indexes of the points in the points vector.
     */
    virtual void setConnectivity(std::vector< std::vector<int> > &
        connectivity_candidate)=0;
    /*! \brief A method for geting all the nodes of the mesh.
     */
    virtual std::vector< std::vector<double> > * getPoints() const =0;
    /*! \brief A method for geting the connectivity of the mesh.
     */
    virtual std::vector< std::vector<int> > * getConnectivity() const =0;
    /*! \brief A method for accessing a specific node of the mesh.
     *
     * @param i index of the point one wishes to acquire. 
     */
    virtual std::vector<double> getPoint(int i) const =0;
    /*! \brief A method for accessing a specific connectivity.
     *
     * @param k index of the connectivity one wishes to acquire. 
     */
    virtual void populatePolytope(int k, Polytope * result) const =0;
    /*! \brief A method for accessing the connectivity of a node.
     *
     * @param i index of the point of which one wishes to get the connectivity 
     * information. 
     */
    virtual std::vector<int> getPointConnectivity(int i) const =0;
    // Helper functions
    virtual int getNumberPoints() const =0;
    virtual int getNumberUnits() const =0;
    virtual int getDimension() const =0;
  protected:
    /*!
     * A pointer to the vector of coordinates describing the nodes of the mesh.
     */
    std::vector< std::vector<double> > points;
    /*!
     * A pointer to the vector of vectors of indexes of points holding the 
     * connectivity information of the mesh.
     */
    std::vector< std::vector<int> > connectivity;
    /*!
     * The dimension of the space.
     */
    int dimension;
};

} //hfox

#endif
