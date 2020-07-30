#include "NabUU.h"

namespace hfox{

void NabUU::setSolution(const std::vector<EVector> & sol){
  solNodes = sol;
  int spaceDim = sol[0].size();
  sols.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  EMatrix solMat(spaceDim, refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solMat.col(i) = sol[i];
  }
  const std::vector<double> * pshape;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    pshape = &(ipShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    sols[i] = solMat*shapes;
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
  EVector locMeasure(refEl->getDimension());
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  EVector gradPhi(nDOFsPerNode);
  EVector testVec(nDOFsPerNode);
  EMatrix buffMat(nDOFsPerNode, nDOFsPerNode);
  EMatrix derivTensor(nDOFsPerNode, nDOFsPerNode);
  double val;
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
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
            op(l*nDOFsPerNode + j, k*nDOFsPerNode + i) -= buffMat.row(i)*derivTensor.col(i);
          }
        }
      }
    }
  }
  for(int i = 0; i < nNodesEl; i++){
    rhs.segment(i*nDOFsPerNode, (i+1)*nDOFsPerNode) = solNodes[i]/2.0;
  }
  rhs = op*rhs;
};//assemble

};//hfox
