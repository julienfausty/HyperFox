#include "Diffusion.h"

namespace hfox{

void Diffusion::assemble(const std::vector< double > & dV, 
    const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Diffusion", "assemble", "cannot assemble before allocating."));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  const std::vector< std::vector< std::vector<double> > > * ipDerivShape;
  ipDerivShape = refEl->getIPDerivShapeFunctions();
  const std::vector<double> * ipWeights = refEl->getIPWeights();
  int elDim = refEl->getDimension();
  int nNodes = refEl->getNumNodes();
  // dV(invJ^T invJ):(\partial\phi_i (\partial \phi_j)^T))(ip)
  std::vector<EVector> invJTinvJs(refEl->getNumIPs(), EVector(elDim*(elDim+1)/2));
  if(diffTensor.size() == 0){
    std::transform(invJacobians.begin(), invJacobians.end(), invJTinvJs.begin(), 
        [this](EMatrix invJ){return symMat2Vec(invJ.transpose() * invJ);});
  } else {
    for(int i = 0; i < refEl->getNumIPs(); i++){
      if(diffTensor[i].cols() == 1){
        invJTinvJs[i] = symMat2Vec(invJacobians[i].transpose() * diffTensor[i](0,0) * invJacobians[i]);
      } else{
        invJTinvJs[i] = symMat2Vec(invJacobians[i].transpose() * diffTensor[i].transpose() * invJacobians[i]);
      }
    }
  }
  EMatrix pPhipPhiT(elDim*(elDim+1)/2, nNodes);
  EMatrix bufferMat(nNodes, nNodes);
  for(int i = 0; i < refEl->getNumIPs(); i++){
    for(int j = 0; j < nNodes; j++){
      Eigen::Map<const EVector> pPhi_j((*ipDerivShape)[i][j].data(), (*ipDerivShape)[i][j].size());
      for(int k = 0; k < nNodes; k++){
        Eigen::Map<const EVector> pPhi_k((*ipDerivShape)[i][k].data(), (*ipDerivShape)[i][k].size());
        bufferMat = pPhi_k * pPhi_j.transpose();
        pPhipPhiT.col(k) = (symMat2Vec(bufferMat + bufferMat.transpose()))/2.0;
      }
      op.block(j, 0, 1, nNodes) += dV[i]*(invJTinvJs[i].transpose() * pPhipPhiT);
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
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

void Diffusion::setDiffTensor(const std::vector<EMatrix> & diffCoeff){
  diffTensor.resize(refEl->getNumIPs(), EMatrix::Zero(diffCoeff[0].rows(), diffCoeff[0].cols()));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  for(int i = 0 ; i < refEl->getNumIPs(); i++){
    for(int j = 0; j < refEl->getNumNodes(); j++){
      diffTensor[i] += (*ipShapes)[i][j]*diffCoeff[j];
    }
  }
};//setDiffTensor

}//hfox
