#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <vector>
#include <algorithm>
#include <mpi.h>
#include "Mesh.h"
#include "Field.h"
#include "ErrorHandle.h"
#include "Modifier.h"
#include "Utils.h"

namespace hfox{

/*!
 * \brief interface class that deals with partitioning the mesh and fields
 *
 * An interface for partitioning problems in HyperFox. Can design interfaces to external libraries from here.
 */

class Partitioner{

  public:
    /*!
     * \brief virtual destructor
     */
    virtual ~Partitioner();
    /*!
     * \brief method for initializing the inner variables
     */
    virtual void initialize();
    /*!
     * \brief method for computing the partition, populates the impotMap and exportMap with relevant partition information
     */
    virtual void computePartition()=0;
    /*!
     * \brief method that updates both the mesh and fields with current computed partition using the importMap/exportMap computed by the computePartition method
     */
    virtual void update();
    /*!
     * \brief method that updates the shared information on each partition
     */
    virtual void updateSharedInformation();
    /*!
     * \brief set the mesh pointer
     * @param pmesh pointer to the mesh
     */
    virtual void setMesh(Mesh * pmesh);
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
    virtual int getRank() const {return rank;};
    /*!
     * \brief get the total number of nodes
     */
    virtual int getTotalNumberNodes() const;
    /*!
     * \brief get the total number of cells
     */
    virtual int getTotalNumberEls() const;
    /*!
     * \brief get the total number of faces
     */
    virtual int getTotalNumberFaces() const;
    /*!
     * \brief read only access to the sharedFaceList
     */
    virtual const std::vector<int> * getSharedFaceList() const{return &sharedFaceList;};
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
     * \brief get the global indeces of a list of nodes from local indeces
     * @param loc list of local indeces of nodes
     * @param glob list to fill with global indeces
     */
    virtual void local2GlobalNodeSlice(const std::vector<int> & loc, std::vector<int> * glob) const;
    /*!
     * \brief get the global indeces of a list of faces from local indeces
     * @param loc list of local indeces of faces
     * @param glob list to fill with global indeces
     */
    virtual void local2GlobalFaceSlice(const std::vector<int> & loc, std::vector<int> * glob) const;
    /*!
     * \brief get the global indeces of a list of elements from local indeces
     * @param loc list of local indeces of elements
     * @param glob list to fill with global indeces
     */
    virtual void local2GlobalElementSlice(const std::vector<int> & loc, std::vector<int> * glob) const;
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
    /*!
     * \brief get the local indeces of a slice of nodes from global indeces
     * @param glob global indeces of nodes
     * @param loc list to fill with local indeces
     */
    virtual void global2LocalNodeSlice(const std::vector<int> & glob, std::vector<int> * loc) const;
    /*!
     * \brief get the local indeces of a slice of faces from global indeces
     * @param glob global indeces of faces
     * @param loc list to fill with local indeces
     */
    virtual void global2LocalFaceSlice(const std::vector<int> & glob, std::vector<int> * loc) const;
    /*!
     * \brief get the local indeces of a slice of elements from global indeces
     * @param glob global indeces of elements
     * @param loc list to fill with local indeces
     */
    virtual void global2LocalElementSlice(const std::vector<int> & glob, std::vector<int> * loc) const;
    /*!
     * \brief get pointer to the node id list
     */
    virtual const std::vector<int> * getNodeIds() const{return &nodeIDs;};
    /*!
     * \brief get pointer to the face id list
     */
    virtual const std::vector<int> * getFaceIds() const{return &faceIDs;};
    /*!
     * \brief get pointer to the cell id list
     */
    virtual const std::vector<int> * getCellIds() const{return &elementIDs;};
  protected:
    /*!
     * \brief helper function for multiplying indexes
     * @param dim number of objects per idex
     * @param idexes the list of indexes
     * @param multiIndexes the list of multiplied indexes
     */
    static void multiplyIndexes(int dim, const std::vector<int> * indexes, std::vector<int> * multiIndexes);
    /*!
     * \brief a method for computing the shared faces of the partition
     */
    void computeSharedFaces();
    /*!
     * \brief a method for empyting out the import and export maps
     */
    void emptyImportExportMaps();
    /*!
     * \brief transmit mesh and field data based on what is indicated in the import and export maps
     * @param iRecvBuffer a map from process to integer data to recieve
     * @param dRecvBuffer a map from process to double data to recieve
     */
    void transmitData(std::map<int, std::vector<int> > * iRecvBuffer, std::map<int, std::vector<double> > * dRecvBuffer, std::map<FieldType, std::vector<Field*> > & fieldMap);
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
     * \brief a list of shared faces of the partition with structure [[globIdFace, rankOtherPartition, globIdAdjacentEl]]
     */
    std::vector<int> sharedFaceList;
    /*!
     * \brief number of partitions
     */
    int nPartitions = 1;
    /*!
     * \brief the index of this partition
     */
    int rank = 0;
    /*!
     * \brief boolean tracking the intialization of the class
     */
    bool initialized = 0;
    /*!
     * \brief a pointer to a map holding the indexes of the cells, faces and nodes to export to another partition exportMap[procToExport][TypeOfObj][index]
     */
    std::map<int, std::map< FieldType, std::vector<int> > > exportMap;
    /*!
     * \brief a pointer to a map holding the indexes of the cells, faces and nodes to import from another partition importMap[procToImport][TypeOfObj][index]
     */
    std::map<int, std::map< FieldType, std::vector<int> > > importMap;

};//Partitioner

}//hfox

#endif//PARTITIONER_H
