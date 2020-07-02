#ifndef ZOLTANOPTS_H
#define ZOLTANOPTS_H

#include <string>
#include<zoltan_cpp.h>

namespace hfox{

/*!
 * \brief Options structure for the ZoltanPartitioner
 */

struct ZoltanOpts{
  /*!
   * \brief the lb_method one wishes to configure Zoltan with
   */
  std::string method = "GRAPH";
  /*!
   * \brief the hypergraph package one wishes to configure Zoltan with
   */
  std::string graphPackage = "PHG";
  /*!
   * \brief the approach one wishes to configure Zoltan with
   */
  std::string approach = "PARTITION";
  /*!
   * \brief configures the debug/verbosity behavior of Zoltan
   */
  std::string debugLevel = "1";
  /*!
   * \brief argument passed to Zoltan_Initialize
   */
  int argc = 0;
  /*!
   * \brief argument passed to Zoltan_Initialize
   */
  char ** argv = NULL;
};


/*!
 * \brief A structure for returning the results from Zoltan's LB_Partition method
 */
struct ZoltanChanges{
  int changes;
  int num_gid_entries;
  int num_lid_entries;
  int num_import;
  ZOLTAN_ID_PTR import_global_ids;
  ZOLTAN_ID_PTR import_local_ids;
  int * import_procs;
  int * import_to_part;
  int num_export;
  ZOLTAN_ID_PTR export_global_ids;
  ZOLTAN_ID_PTR export_local_ids;
  int * export_procs;
  int * export_to_part;
};

};

#endif//ZOLTANOPTS_H
