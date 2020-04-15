#include "DirichletCondition.h"

namespace hfox{

void DirichletCondition::assemble(const std::vector< double > & detJacobians, const std::vector< EMatrix > & invJacobians){
  op = EMatrix::Identity(op.rows(), op.cols());
};//assemble

}//hfox
