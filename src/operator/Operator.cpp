#include "Operator.h"

namespace hfox{

void Operator::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser < 1){
    throw(ErrorHandle("Operator", "allocate", "the number of DOFs per node must be at least one"));
  }
  nDOFsPerNode = nDOFsPerNodeUser;
  op = Eigen::MatrixXd::Zero(nDOFsPerNode * refEl->getNumNodes(), nDOFsPerNode * refEl->getNumNodes());
  allocated = 1;
};//allocate

std::vector<EMatrix> Operator::calcJacobians(const std::vector< std::vector<double> > & points,
    const ReferenceElement * referenceEl){
  if(points.size() != referenceEl->getNumNodes()){
    throw(ErrorHandle("Operator", "calcJacobians", 
          "the number of nodes passed to the method must be equal to the " 
          "number of nodes in the reference element."));
  }
  int dimMeshSpace = points[0].size();
  int dimRefSpace = referenceEl->getDimension();
  if(dimMeshSpace < dimRefSpace){
    throw(ErrorHandle("Operator", "calcJacobians", 
          "the dimension of the mesh space cannot be smaller than the "
          "topological dimension of the reference element."));
  }
  std::vector<EMatrix> jacobians(referenceEl->getNumIPs(), Eigen::MatrixXd::Zero(dimRefSpace, dimMeshSpace));
  const std::vector< std::vector< std::vector<double> > > * ipDerivShape;
  ipDerivShape = referenceEl->getIPDerivShapeFunctions();
  for(int i = 0; i < referenceEl->getNumNodes(); i++){
    Eigen::Map<const EMatrix> elNodes(points[i].data(), 1, dimMeshSpace);
    for(int j = 0; j < referenceEl->getNumIPs(); j++){
      Eigen::Map<const EVector> derivs((*ipDerivShape)[j][i].data(), dimRefSpace);
      jacobians[j] += derivs * elNodes;
    }
  }
  return jacobians;
};//calcJacobians

std::vector<EMatrix> Operator::calcInvJacobians(const std::vector< EMatrix > & jacobians){
  int nJacs = jacobians.size();
  if(nJacs != 0){
    int rowsJacs = jacobians[0].rows();
    int colsJacs = jacobians[0].cols();
    std::vector<EMatrix> invJacs(nJacs, EMatrix::Zero(rowsJacs, colsJacs));
    if(rowsJacs == colsJacs){
      std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){return jac.inverse();});
    } else {
      std::transform(jacobians.begin(), jacobians.end(), invJacs.begin(), [](EMatrix jac){
          return jac.completeOrthogonalDecomposition().pseudoInverse();});
    }
    return invJacs;
  } else {
    return std::vector<EMatrix>(0);
  }
};//calcInvJacobians

std::vector<double> Operator::calcDetJacobians(const std::vector< EMatrix > & jacobians){
  int nJacs = jacobians.size();
  if(nJacs != 0){
    int rowsJacs = jacobians[0].rows();
    int colsJacs = jacobians[0].cols();
    std::vector<double> detJacs(nJacs, 0.0);
    if(rowsJacs == colsJacs){
      std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){return jac.determinant();});
    } else {
      std::transform(jacobians.begin(), jacobians.end(), detJacs.begin(), [](EMatrix jac){
          return std::sqrt((jac*jac.transpose()).determinant());});
    }
    return detJacs;
  } else {
    return std::vector<double>(0);
  }
};//calcInvJacobians

}//hfox
