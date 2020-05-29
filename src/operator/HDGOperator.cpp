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

void HDGOperator::multiplyDOFs(){
  int nNodes = op.cols()/nDOFsPerNode;
  EMatrix buff = op.block(0, 0, nNodes, nNodes);
  op.block(0,0,nNodes,nNodes) = EMatrix::Zero(nNodes, nNodes);
  for(int i = 0; i < nNodes; i++){
    for(int j = 0; j < nNodes; j++){
      for(int k = 0; k < nDOFsPerNode; k++){
        op(i*nDOFsPerNode + k, j*nDOFsPerNode + k) = buff(i, j);
      }
    }
  }
};//multiplyDOFs

}//hfox
