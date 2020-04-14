#include "Operator.h"

namespace hfox{

void Operator::allocate(int nDOFsPerNodeUser){
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
  std::vector<EMatrix> jacobians(referenceEl->getNumIPs(), Eigen::MatrixXd::Zero(dimMeshSpace, dimRefSpace));
  const std::vector< std::vector< std::vector<double> > > * ipDerivShape;
  ipDerivShape = referenceEl->getIPDerivShapeFunctions();
  for(int i = 0; i < referenceEl->getNumNodes(); i++){
    Eigen::Map<const EVector> elNodes(points[i].data(), dimMeshSpace);
    for(int j = 0; j < referenceEl->getNumIPs(); j++){
      Eigen::Map<const EMatrix> derivs((*ipDerivShape)[j][i].data(), 1, dimRefSpace);
      jacobians[j] += elNodes * derivs;
    }
  }
  return jacobians;
};

}//hfox
