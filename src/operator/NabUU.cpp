#include "NabUU.h"

namespace hfox{

void NabUU::setSolution(const std::vector<EVector> & sol){
  int spaceDim = sol[0].size();
  solNodes.resize(refEl->getNumNodes(), EVector::Zero(spaceDim));
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solNodes[i] = sol[i];
  }
  sols.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  gradSols.resize(refEl->getNumIPs(), EMatrix::Zero(spaceDim, spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  EMatrix solMat(spaceDim, refEl->getNumNodes());
  EMatrix gradMat(refEl->getNumNodes(), spaceDim);
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solMat.col(i) = sol[i];
  }
  const std::vector<double> * pshape;
  const std::vector< std::vector<double> > * derivShape;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    pshape = &(ipShapes->at(i));
    derivShape = &(ipDerivShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    sols[i] = solMat*shapes;
    for(int k = 0; k < refEl->getNumNodes(); k++){
      gradMat.row(k) = EMap<const EVector>(derivShape->at(k).data(), derivShape->at(k).size()).transpose();
    }
    gradSols[i] = solMat * gradMat;
  }
};//setSolution

void NabUU::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("NabUU", "assemble", "the operator should be allocated before being assembled."));
  }
  if(sols.size() == 0){
    throw(ErrorHandle("NabUU", "assemble", "the solution must be set before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  rhs = EVector::Zero(rhs.size());
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  EVector gradPhi(nDOFsPerNode);
  EVector testVec(nDOFsPerNode);
  EMatrix buffMat(nDOFsPerNode, nDOFsPerNode);
  EMatrix derivTensor(nDOFsPerNode, nDOFsPerNode);
  EMatrix gradMat;
  double divSol;
  double val;
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    gradMat = gradSols[ip] * invJacobians[ip].transpose() * dV[ip];
    divSol = gradMat.trace();
    for(int k = 0; k < nNodesEl; k++){
      for(int i = 0; i < nDOFsPerNode; i++){
        testVec = EVector::Zero(nDOFsPerNode);
        testVec[i] = dV[ip] * shapes->at(k);
        buffMat = (testVec) * sols[ip].transpose();
        for(int l = 0; l < nNodesEl; l++){
          gradPhi = invJacobians[ip]*EMap<const EVector>(derivShapes->at(l).data(), derivShapes->at(l).size());
          for(int j = 0; j < nDOFsPerNode; j++){
            derivTensor = EMatrix::Zero(nDOFsPerNode, nDOFsPerNode);
            derivTensor.col(j) = gradPhi;
            derivTensor += derivTensor.transpose();
            op(k*nDOFsPerNode + i, l*nDOFsPerNode + j) += buffMat.row(i)*derivTensor.col(i) + gradMat(j, i)*shapes->at(k)*shapes->at(l);
            if(i == j){
              op(k*nDOFsPerNode + i, l*nDOFsPerNode + j) += divSol * shapes->at(k) * shapes->at(l);
            }
          }
        }
      }
    }
  }
  for(int i = 0; i < nNodesEl; i++){
    for(int j = 0; j < nDOFsPerNode; j++){
      rhs[i*nDOFsPerNode + j] = solNodes[i][j]/2.0;
    }
  }
  rhs = op*rhs;
};//assemble

};//hfox
