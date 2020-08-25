#include "HDGNonLinearOperator.h"

namespace hfox{


void HDGNonLinearOperator::allocate(int nDOFsPerNodeUser){
  HDGOperator::allocate(nDOFsPerNodeUser);
  rhs = EVector::Zero(op.cols());
};//allocate

}
