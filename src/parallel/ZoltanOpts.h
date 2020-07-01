#ifndef ZOLTANOPTS_H
#define ZOLTANOPTS_H

#include <string>

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
   * \brief argument passed to Zoltan_Initialize
   */
  int argc = 0;
  /*!
   * \brief argument passed to Zoltan_Initialize
   */
  char ** argv = NULL;
};

};

#endif//ZOLTANOPTS_H
