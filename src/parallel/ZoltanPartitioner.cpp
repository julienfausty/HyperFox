#include "ZoltanPartitioner.h"

namespace hfox{

void ZoltanPartitioner::initialize(){
  Partitioner::initialize();
  float version;
  Zoltan_Initialize(myOpts.argc, myOpts.argv, &version);
};//initialize

};//hfox
