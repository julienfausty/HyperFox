#include "Convection.h"

namespace hfox{

void Convection::setVelocity(const std::vector<EVector> & velocity){
  int spaceDim = velocity[0].size();
  vels.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  EMatrix velMat(spaceDim, refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    velMat.col(i) = velocity[i];
  }
  for(int i = 0; i < refEl->getNumIPs(); i++){
    const std::vector<double> * pshape = &(ipShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    vels[i] = velMat*shapes;
  }
};//setVelocity

void Convection::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Convection", "assemble", "the operator should be allocated before being assembled."));
  }
  if(vels.size() == 0){
    throw(ErrorHandle("Convection", "assemble", "the velocity must be set before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  EVector locMeasure(refEl->getDimension());
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  const std::vector<double> * partialPhi;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    locMeasure = (((vels[i]).transpose() * invJacobians[i])*(dV[i])).transpose();
    shapes = &(ipShapes->at(i));
    derivShapes = &(ipDerivShapes->at(i));
    for(int k = 0; k < refEl->getNumNodes(); k++){
      for(int l = 0; l < refEl->getNumNodes(); l++){
        partialPhi = &(derivShapes->at(l));
        op(k, l) += locMeasure.dot(Eigen::Map<const EVector>(partialPhi->data(), partialPhi->size())) * (shapes->at(k));
      }
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
};//assemble

}//hfox
