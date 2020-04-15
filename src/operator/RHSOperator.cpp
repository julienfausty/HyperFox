#include "RHSOperator.h"

namespace hfox{

void RHSOperator::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser < 1){
    throw(ErrorHandle("RHSOperator", "allocate", "the number of DOFs per node must be at least one"));
  }
  nDOFsPerNode = nDOFsPerNodeUser;
  op = Eigen::MatrixXd::Zero(nDOFsPerNode * refEl->getNumNodes(), 1);
  allocated = 1;
};//allocate

}//hfox
