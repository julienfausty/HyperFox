#include "Diffusion.h"

namespace hfox{

void Diffusion::assemble(const std::vector< double > & detJacobians, 
    const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Diffusion", "assemble", "cannot assemble before allocating."));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  // w(ip)(det(J)(invJ^T invJ):(\partial\phi_i (\partial \phi_j)^T))(ip)
  std::vector<EVector> invJTinvJs(refEl->getNumIPs());
  std::transform(invJacobians.begin(), invJacobians.end(), invJTinvJs.begin(), 
      [this](EMatrix invJ){return symMat2Vec(invJ.transpose() * invJ);});
  const std::vector< std::vector< std::vector<double> > > * ipDerivShape;
  ipDerivShape = refEl->getIPDerivShapeFunctions();
  const std::vector<double> * ipWeights = refEl->getIPWeights();
  int elDim = refEl->getDimension();
  int nNodes = refEl->getNumNodes();
  for(int i = 0; i < refEl->getNumIPs(); i++){
    EMatrix pPhipPhiT(elDim*(elDim+1)/2, nNodes);
    for(int j = 0; j < nNodes; j++){
      Eigen::Map<const EVector> pPhi_j((*ipDerivShape)[i][j].data(), (*ipDerivShape)[i][j].size());
      for(int k = 0; k < nNodes; k++){
        Eigen::Map<const EVector> pPhi_k((*ipDerivShape)[i][k].data(), (*ipDerivShape)[i][k].size());
        pPhipPhiT.col(k) = symMat2Vec(pPhi_k * pPhi_j.transpose());
      }
      op.block(j, 0, 1, nNodes) += ((*ipWeights)[i])*detJacobians[i]*(invJTinvJs[i].transpose() * pPhipPhiT);
    }
  }
  if(nDOFsPerNode > 1){
    EMatrix buff = op.block(0, 0, nNodes, nNodes);
    op.block(0,0,nNodes,nNodes) = EMatrix::Zero(nNodes, nNodes);
    for(int i = 0; i < nNodes; i++){
      for(int j = 0; j < nNodes; j++){
        for(int k = 0; k < nDOFsPerNode; k++){
          op(i*nDOFsPerNode + k, j*nDOFsPerNode + k) = buff(i, j);
        }
      }
    }
  }
};//assemble

EVector Diffusion::symMat2Vec(const EMatrix & symMat){
  int n = symMat.cols();
  EVector res(n*(n+1)/2);
  double sqrt2 = std::pow(2.0, 1.0/2.0);
  int index = 0;
  for(int i = 0; i < n ; i++){
    for(int j = i; j < n; j++){
      res[index] = symMat(i, j);
      if(i != j){
        res[index] *= sqrt2;
      }
      index += 1;
    }
  }
  return res;
};//symMat2Vec

}//hfox
