#include "HDGOperator.h"

namespace hfox{

void HDGOperator::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser < 1){
    throw(ErrorHandle("HDGOperator", "allocate", "the number of DOFs per node must be at least one"));
  }
  nDOFsPerNode = nDOFsPerNodeUser;
  //(Solution (nNodes) + Flux (nNodes * dim) + Trace (nNodesFace * nFaces))
  int n = nDOFsPerNode * (refEl->getNumNodes()*(refEl->getDimension() + 1) 
      + (refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces()));
  op = EMatrix::Zero(n, n);
  allocated = 1;
};//allocate

}//hfox
