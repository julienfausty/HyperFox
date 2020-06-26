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
    ZoltanPartitioner(){;};
    /*!
     * \brief option constructor
     */
    ZoltanPartitioner(ZoltanOpts opts){setOptions(opts);};
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
     *
     * @param weights parameters controlling the different weights given to each partition
     */
    virtual void computePartition(std::vector<double> weights = std::vector<double>());
    /*!
     * \brief method that updates both the mesh and fields with current computed partition
     */
    virtual void update();

  private:
    /*!
     * \brief the special options for the Partitioner
     */
    ZoltanOpts myOpts;

};//ZoltanPartitioner

}//hfox

#endif//ZOLTANPARTITIONER_H
