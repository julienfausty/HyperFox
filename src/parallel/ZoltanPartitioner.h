#ifndef ZOLTANPARTITIONER_H
#define ZOLTANPARTITIONER_H

#include <zoltan_cpp.h>
#include "Partitioner.h"
#include "ZoltanOpts.h"

namespace hfox{

/*!
 * \brief interface with zoltan partitioner
 */

class ZoltanPartitioner : public Partitioner{

  public:
    /*!
     * \brief Constructor
     */
    ZoltanPartitioner(Mesh * pMesh){setMesh(pMesh);};
    /*!
     * \brief destructor
     */
    ~ZoltanPartitioner();
    /*!
     * \brief option constructor
     */
    ZoltanPartitioner(Mesh * pMesh, ZoltanOpts opts){setMesh(pMesh), setOptions(opts);};
    /*!
     * \brief set the options for the partitioner
     */
    void setOptions(ZoltanOpts opts){myOpts = opts;};
    /*!
     * \brief method for initializing the inner variables
     */
    virtual void initialize();
    /*!
     * \brief method for computing the partition
     */
    virtual void computePartition();
    /*!
     * \brief method that updates both the mesh and fields with current computed partition
     */
    virtual void update();

  private:
    /*!
     * \brief the special options for the Partitioner
     */
    ZoltanOpts myOpts;
    /*!
     * \brief the Zoltan object from the Zoltan library
     */
    Zoltan * zObj = NULL;
    /*!
     * \brief a struct for storing the results of the changes to the partition
     */
    ZoltanChanges zChanges;
    /*!
     * \brief static function for getting the number of cells on partition (for Zoltan)
     *
     * @param data this is supposed to be the mesh
     * @param ierr an arror code for Zoltan
     */
    static int getNumberCells(void * data, int * ierr);
    /*!
     * \brief static function for getting the list of cells on partition (for Zoltan)
     *
     * @param data this is the list of element ids
     * @param num_gid_entries the dimension of a global index (should always be one here)
     * @param num_lid_entries the dimension of a local index (should always be one here)
     * @param global_ids the return array for the global indexes of the cells on the partition
     * @param local_ids the return array for the local indexes of the cells on the partition
     * @param wgt_dim the dimension of the weight object
     * @param obj_weights a return array with weights given to all the cells (not supported here yet)
     * @param ierr an arror code for Zoltan
     */
    static void getListCells(void * data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float * obj_weights, int * ierr);
    /*!
     * \brief static function for getting the number of neighbors for a given list of global ids of cells
     *
     * @param data this is a std::pair object with pointers to the parallel mesh and the partitioner
     * @param num_gid_entries the dimension of a global index (should always be one here)
     * @param num_lid_entries the dimension of a local index (should always be one here)
     * @param num_obj number of objects in the list of cells
     * @param global_ids the array of global ids of cells that one should use
     * @param local_ids the array of local ids of cells that one should look for
     * @param num_edges the return array of number of edges for each index
     * @param ierr an arror code for Zoltan
     */
    static void getNumNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, int * ierr);
    /*!
     * \brief static function for getting the global ids of the neighbors of a given list of global ids of cells
     *
     * @param data this is a std::pair object with pointers to the parallel mesh and the partitioner
     * @param num_gid_entries the dimension of a global index (should always be one here)
     * @param num_lid_entries the dimension of a local index (should always be one here)
     * @param num_obj number of objects in the list of cells
     * @param global_ids the array of global ids of cells that one should use
     * @param local_ids the array of local ids of cells that one should look for
     * @param num_edges the number of neighbors for each cell
     * @param nbor_global_ids the return array to fill with neighbors for each cell
     * @param wgt_dim the dimension of the weight object
     * @param ewgts a return array with weights given to all the edges (not supported here yet)
     * @param ierr an arror code for Zoltan
     */
    static void getListNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, ZOLTAN_ID_PTR nbor_global_ids, int * nbor_procs, int wgt_dim, float * ewgts, int * ierr);


};//ZoltanPartitioner

}//hfox

#endif//ZOLTANPARTITIONER_H
