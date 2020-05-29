#include "HDGConvection.h"

namespace hfox{

HDGConvection::HDGConvection(const ReferenceElement * re) : HDGOperator(re){
  mass = new Mass(re);
};//constructor

HDGConvection::~HDGConvection(){
  delete mass;
};//destructor

void HDGConvection::allocate(int nDOFsPerNodeUser){
  mass->allocate(1);
  HDGOperator::allocate(nDOFsPerNodeUser);
};//allocate

void HDGConvection::setVelocity(const std::vector<EVector> & vels){
  int spaceDim = vels[0].size();
  velocities.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  EMatrix velMat(spaceDim, refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    velMat.col(i) = vels[i];
  }
  const std::vector<double> * pshape;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    pshape = &(ipShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    velocities[i] = velMat*shapes;
  }
};//setVelocity

void HDGConvection::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGConvection", "assemble", "the operator should be allocated before being assembled."));
  }
  if(velocities.size() == 0){
    throw(ErrorHandle("HDGConvection", "assemble", "the velocity must be set before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  int dimSpace = refEl->getDimension();
  int nNodes = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  std::vector<double> locMeasure(nIPs, 0.0);
  for(int i = 0; i < dimSpace; i++){
    for(int k = 0; k < nIPs; k++){
      locMeasure[k] = dV[k]*velocities[k][i];
    }
    mass->assemble(locMeasure, invJacobians);
    for(int k = 0; k < nNodes; k++){
      op.block(0, nNodes + (dimSpace)*k + i, nNodes, 1) = mass->getMatrix()->col(k);
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
};//assemble

}//hfox
