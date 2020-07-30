#include "NonLinearOperator.h"

namespace hfox{

void NonLinearOperator::allocate(int nDOFsPerNodeUser){
  Operator::allocate(nDOFsPerNodeUser);
  rhs = EVector::Zero(nDOFsPerNode * (refEl->getNumNodes()));
};//allocate

};//hfox
