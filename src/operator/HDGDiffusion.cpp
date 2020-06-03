#include "HDGDiffusion.h"

namespace hfox{

HDGDiffusion::HDGDiffusion(const ReferenceElement * re) : HDGOperator(re){
  convection = new Convection(re);
  faceMass = new Mass(re->getFaceElement());
};//constructor

HDGDiffusion::~HDGDiffusion(){
  delete convection;
  delete faceMass;
};//destructor

void HDGDiffusion::allocate(int nDOFsPerNodeUser){
  convection->allocate(1);
  faceMass->allocate(1);
  HDGOperator::allocate(nDOFsPerNodeUser);
};//allocate

void HDGDiffusion::setFromBase(const std::vector<EVector> * ns){
  if(ns->at(0).size() != refEl->getDimension()){
    throw(ErrorHandle("HDGDiffusion", "setFromBase", "the dimension of the normals must be equal to the dimension of the reference space"));
  }
  if(ns->size() != refEl->getNumFaces() * (refEl->getFaceElement()->getNumIPs())){
    throw(ErrorHandle("HDGDiffusion", "setFromBase", "the number of normals should be the number of faces times the number of integration points per face"));
  }
  normals = ns;
};

void HDGDiffusion::setDiffusionTensor(const std::vector<EMatrix> & diffTensor){
  int nNodes = refEl->getNumNodes();
  if(diffTensor.size() != nNodes){
    throw(ErrorHandle("HDGDiffusion", "setDiffusionTensor", "the number of matrices should be the number nodes in the element"));
  }
  int dim = refEl->getDimension();
  if((diffTensor[0].rows() != dim) or (diffTensor[0].cols() != dim)){
    throw(ErrorHandle("HDGDiffusion", "setDiffusionTensor", "the diffusion tensor matrices should be square and have a size of the reference dimension times the number of DOFs per node"));
  }
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nIPs = fEl->getNumIPs();
  int nNodesPFc = fEl->getNumNodes();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * ipShapes = fEl->getIPShapeFunctions();
  Ds.resize(nFaces*nIPs, EMatrix::Zero(dim, dim));
  EMatrix nodeDs(nNodesPFc, dim*dim);
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int i = 0; i < nNodesPFc; i++){
      nodeDs.row(i) = EMap<const EMatrix>(diffTensor[faceNodeMap->at(iFace)[i]].data(), 1, dim*dim);
    }
    for(int i = 0; i < nIPs; i++){
      EMap<const EVector> shape(ipShapes->at(i).data(), ipShapes->at(i).size());
      EMap<EMatrix> res(Ds[iFace*nIPs + i].data(), 1, dim*dim);
      res = shape.transpose()*nodeDs;
    }
  }
  myDiffTensor = diffTensor;
};//setDiffusionTensor

void HDGDiffusion::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGDiffusion", "assemble", "must allocated the operator before assembling"));
  }
  if(normals == NULL){
    throw(ErrorHandle("HDGDiffusion", "assemble", "must set the normals from the base before assembling"));
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
  if(Ds.size() == 0){
    Ds.resize(nFaces * nIPsPFc, EMatrix::Identity(dim, dim));
    myDiffTensor.resize(nNodesEl, EMatrix::Identity(dim, dim));
  }
  std::vector<EVector> velocity(nNodesEl, EVector::Zero(dim));
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < nNodesEl; j++){
      velocity[j] = myDiffTensor[j].col(i);
    }
    convection->setVelocity(velocity);
    convection->assemble(dV, invJacobians);
    for(int j = 0; j < nNodesEl; j++){
      op.block(0, nNodesEl + j*dim + i, nNodesEl, 1) += (convection->getMatrix()->row(j)).transpose();
    }
  }
  std::vector<EMatrix> faceInvJacs(nIPsPFc);
  std::vector<double> locMeasure(nIPsPFc, 0.0);
  std::vector<EVector> diffNorms(nIPsPFc, EVector::Zero(dim));
  int offset = 0, facesOff = 0;
  for(int i = 0; i < nFaces; i++){
    facesOff = i*nIPsPFc;
    offset = nIPsEl + facesOff;
    std::copy(invJacobians.begin() + offset, invJacobians.begin() + offset + nIPsPFc, faceInvJacs.begin());
    for(int j = 0; j < nIPsPFc; j++){
      diffNorms[j] = Ds[facesOff + j].transpose()*(normals->at(facesOff + j));
    }
    for(int d = 0; d < dim; d++){
      for(int j = 0; j < nIPsPFc; j++){
        locMeasure[j] = dV[offset + j]*(diffNorms[j][d]);
      }
      faceMass->assemble(locMeasure, faceInvJacs);
      for(int j = 0; j < nNodesFace; j++){
        for(int k = 0; k < nNodesFace; k++){
          op(faceNodeMap->at(i)[k], nNodesEl + dim*(faceNodeMap->at(i)[j]) + d) -= (*(faceMass->getMatrix()))(k, j);
          op(nNodesEl*(dim+1) + i*nNodesFace + k, nNodesEl + dim*(faceNodeMap->at(i)[j]) + d) -= (*(faceMass->getMatrix()))(k, j);
        }
      }
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
};//assemble

}//hfox
