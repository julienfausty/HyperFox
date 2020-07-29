#include "HDGLinearizedConvection.h"

namespace hfox{

HDGLinearizedConvection::HDGLinearizedConvection(const ReferenceElement * re) : HDGOperator(re){
  convection = new Convection(re);
  faceMass = new Mass(re->getFaceElement());
};//constructor

HDGLinearizedConvection::~HDGLinearizedConvection(){
  delete convection;
  delete faceMass;
};//destructor

void HDGLinearizedConvection::allocate(int nDOFsPerNodeUser){
  convection->allocate(1);
  faceMass->allocate(1);
  HDGOperator::allocate(nDOFsPerNodeUser);
};//allocate

void HDGLinearizedConvection::setFromBase(const std::vector<EVector> * ns){
  if(ns->at(0).size() != refEl->getDimension()){
    throw(ErrorHandle("HDGLinearizedConvection", "setFromBase", "the dimension of the normals must be equal to the dimension of the reference space"));
  }
  if(ns->size() != refEl->getNumFaces() * (refEl->getFaceElement()->getNumIPs())){
    throw(ErrorHandle("HDGLinearizedConvection", "setFromBase", "the number of normals should be the number of faces times the number of integration points per face"));
  }
  normals = ns;
};//setFromBase

void HDGLinearizedConvection::setVelocity(const std::vector<EVector> & vels){
  int nNodes = refEl->getNumNodes();
  if(vels.size() != nNodes){
    throw(ErrorHandle("HDGLinearizedConvection", "setVelocity", "the number of vectors should be the number nodes in the element"));
  }
  int dim = refEl->getDimension();
  if((vels[0].size() != dim)){
    throw(ErrorHandle("HDGLinearizedConvection", "setVelocity", "the velocities should have the same number of components as the reference element dimension."));
  }
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nIPs = fEl->getNumIPs();
  int nNodesPFc = fEl->getNumNodes();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * ipShapes = fEl->getIPShapeFunctions();
  faceVels.resize(nFaces*nIPs, EVector::Zero(dim));
  EMatrix nodeVels(dim, nNodesPFc);
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int i = 0; i < nNodesPFc; i++){
      nodeVels.col(i) = vels[faceNodeMap->at(iFace)[i]];
    }
    for(int i = 0; i < nIPs; i++){
      EMap<const EVector> shape(ipShapes->at(i).data(), ipShapes->at(i).size());
      faceVels[iFace*nIPs + i] = nodeVels*shape;
    }
  }
  velocities = vels;
};//setVelocity

void HDGLinearizedConvection::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGLinearizedConvection", "assemble", "the operator should be allocated before being assembled."));
  }
  if(velocities.size() == 0){
    throw(ErrorHandle("HDGLinearizedConvection", "assemble", "the velocity must be set before assembling"));
  }
  if(normals == NULL){
    throw(ErrorHandle("HDGLinearizedConvection", "assemble", "must set the normals from the base before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  int dim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFace = fEl->getNumNodes();
  int nIPsPFc = fEl->getNumIPs();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes(); 
  convection->setVelocity(velocities);
  convection->assemble(dV, invJacobians);
  op.block(0, 0, nNodesEl, nNodesEl) -= convection->getMatrix()->transpose();
  std::vector<EMatrix> faceInvJacs(nIPsPFc);
  std::vector<double> locMeasure(nIPsPFc, 0.0);
  int offset = 0, facesOff = 0;
  for(int i = 0; i < nFaces; i++){
    facesOff = i*nIPsPFc;
    offset = nIPsEl + facesOff;
    std::copy(invJacobians.begin() + offset, invJacobians.begin() + offset + nIPsPFc, faceInvJacs.begin());
    for(int j = 0; j < nIPsPFc; j++){
      locMeasure[j] = dV[offset + j]*(faceVels[facesOff + j].dot(normals->at(facesOff + j)));
    }
    faceMass->assemble(locMeasure, faceInvJacs);
    offset = nNodesEl*(dim+1) + i*nNodesFace;
    for(int j = 0; j < nNodesFace; j++){
      for(int k = 0; k < nNodesFace; k++){
        op(faceNodeMap->at(i)[k], offset + j) += (*(faceMass->getMatrix()))(k, j);
        op(offset + k, offset + j) += (*(faceMass->getMatrix()))(k, j);
      }
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
};//assemble

}//hfox
