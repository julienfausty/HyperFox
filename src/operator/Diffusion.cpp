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
  int nIPsEl = refEl->getNumIPs();
  const std::vector< std::vector<double> > * derivShapes;
  EMatrix gradShapes(elDim, nNodes);
  std::vector<EMatrix> diffusionTensorsIPs(nIPsEl, EMatrix::Identity(elDim, elDim));
  if(diffTensor.size() != 0){
    for(int ip = 0; ip < nIPsEl; ip++){
      if(diffTensor[ip].cols() == 1){
        diffusionTensorsIPs[ip] = diffTensor[ip](0, 0)*diffusionTensorsIPs[ip];
      } else {
        diffusionTensorsIPs[ip] = diffTensor[ip];
      }
    }
  }
  for(int ip = 0; ip < nIPsEl; ip++){
    derivShapes = &(ipDerivShape->at(ip));
    for(int iN = 0; iN < nNodes; iN++){
      gradShapes.col(iN) = invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), elDim);
    }
    for(int iN = 0; iN < nNodes; iN++){
      for(int jN = 0; jN < nNodes; jN++){
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) += 
            ((diffusionTensorsIPs[ip]*gradShapes.col(iN)).transpose()*gradShapes.col(jN))(0,0)*dV[ip];
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

void Diffusion::setDiffTensor(const std::vector<EMatrix> & diffCoeff){
  diffTensor.resize(0);
  diffTensor.resize(refEl->getNumIPs(), EMatrix::Zero(diffCoeff[0].rows(), diffCoeff[0].cols()));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  for(int i = 0 ; i < refEl->getNumIPs(); i++){
    for(int j = 0; j < refEl->getNumNodes(); j++){
      diffTensor[i] += (*ipShapes)[i][j]*diffCoeff[j];
    }
  }
};//setDiffTensor

}//hfox
