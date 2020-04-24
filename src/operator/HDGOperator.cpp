#include "HDGOperator.h"

namespace hfox{

void HDGOperator::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser < 1){
    throw(ErrorHandle("HDGOperator", "allocate", "the number of DOFs per node must be at least one"));
  }
  nDOFsPerNode = nDOFsPerNodeUser;
  //(Solution (nNodes) + Flux (nNodes * dim) + Trace (nFaceNodes))
  int n = nDOFsPerNode * (refEl->getNumNodes()*(refEl->getDimension() + 2) - refEl->getInnerNodes()->size());
  op = EMatrix::Zero(n, n);
  allocated = 1;
};//allocate

}//hfox
