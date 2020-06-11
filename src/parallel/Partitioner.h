#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <vector>
#include "Mesh.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief interface class that deals with partitioning the mesh and fields
 *
 * An interface for partitioning problems in HyperFox. Can design interfaces to external libraries from here.
 */

class Partitioner{

  public:
    /*!
     * \brief method for initializing the inner variables
     */
    virtual void initialize()=0;
    /*!
     * \brief method for computing the partition
     *
     * @param weights parameters controlling the different weights given to each partition
     */
    virtual void computePartition(std::vector<double> weights = std::vector<double>())=0;
    /*!
     * \brief method that updates both the mesh and fields with current computed partition
     */
    virtual void update()=0;
    /*!
     * \brief set the mesh pointer
     * @param pmesh pointer to the mesh
     */
    virtual void setMesh(Mesh * pmesh){myMesh = pmesh;};
    /*!
     * \brief set the field list
     * @param fieldList list of pointers to the fields
     */
    virtual void setFields(std::vector<Field*> fieldList){fields = fieldList;};
    /*!
     * \brief get number of partitions
     */
    virtual int getNumPartitions() const {return nPartitions;};
    /*!
     * \brief get the index of this partition
     */
    virtual int getThisPartition() const {return thisPartition;};
    /*!
     * \brief get the global index of node from a local index
     * @param loc local index of node
     */
    virtual int local2GlobalNode(int loc) const {return nodeIDs[loc];};
    /*!
     * \brief get the global index of face from a local index
     * @param loc local index of face
     */
    virtual int local2GlobalFace(int loc) const {return faceIDs[loc];};
    /*!
     * \brief get the global index of element from a local index
     * @param loc local index of element
     */
    virtual int local2GlobalEl(int loc) const {return elementIDs[loc];};
    /*!
     * \brief get the local index of node from a global index
     * @param glob global index of node
     */
    virtual int global2LocalNode(int glob) const;
    /*!
     * \brief get the local index of face from a global index
     * @param glob global index of face
     */
    virtual int global2LocalFace(int glob) const;
    /*!
     * \brief get the local index of element from a global index
     * @param glob global index of element
     */
    virtual int global2LocalElement(int glob) const;
  protected:
    /*!
     * \brief pointer to the partitioned mesh
     */
    Mesh * myMesh = NULL;
    /*!
     * \brief list of all fields attached to the mesh
     */
    std::vector<Field*> fields;
    /*!
     * \brief list of global node IDs
     */
    std::vector<int> nodeIDs;
    /*!
     * \brief list of global face IDs
     */
    std::vector<int> faceIDs;
    /*!
     * \brief list of global element IDs
     */
    std::vector<int> elementIDs;
    /*!
     * \brief number of partitions
     */
    int nPartitions = 1;
    /*!
     * \brief the index of this partition
     */
    int thisPartition = 0;
    /*!
     * \brief boolean tracking the intialization of the class
     */
    bool initialized = 0;

};//Partitioner

}//hfox

#endif//PARTITIONER_H
