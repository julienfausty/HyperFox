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
  int nIPsEl = refEl->getNumIPs();
  int nIPsFc = fEl->getNumIPs();
  int nNodesEl = refEl->getNumNodes();
  int nNodesPFc = fEl->getNumNodes();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  Ds.resize(0);
  Ds.resize(nIPsEl + nFaces*nIPsFc, EMatrix::Zero(dim, dim));
  EMatrix nodeDs(nNodesEl, dim*dim);
  for(int i = 0; i < nNodesEl; i++){
    nodeDs.row(i) = EMap<const EMatrix>(diffTensor[i].data(), 1, dim*dim);
  }
  for(int i = 0; i < nIPsEl; i++){
    EMap<const EVector> shape(ipShapes->at(i).data(), ipShapes->at(i).size());
    EMap<EMatrix> res(Ds[i].data(), 1, dim*dim);
    res = shape.transpose()*nodeDs;
  }
  nodeDs.resize(nNodesPFc, dim*dim);
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int i = 0; i < nNodesPFc; i++){
      nodeDs.row(i) = EMap<const EMatrix>(diffTensor[faceNodeMap->at(iFace)[i]].data(), 1, dim*dim);
    }
    for(int i = 0; i < nIPsFc; i++){
      EMap<const EVector> shape(fipShapes->at(i).data(), fipShapes->at(i).size());
      EMap<EMatrix> res(Ds[nIPsEl + iFace*nIPsFc + i].data(), 1, dim*dim);
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
  int spaceDim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesPFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  const std::vector<int> * faceNodes;
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  int lenU = nNodesEl * nDOFsPerNode;
  int lenQ = spaceDim*lenU;
  int lenL = nNodesPFc*nFaces*nDOFsPerNode;
  int offset;
  EVector buffVec(spaceDim);
  EVector buffVec2(spaceDim);
  EVector buffVec3(spaceDim);
  EMatrix buffMat(spaceDim, nDOFsPerNode);
  if(Ds.size() == 0){
    Ds.resize(nIPsEl + nFaces * nIPsFc, EMatrix::Identity(spaceDim, spaceDim));
    myDiffTensor.resize(nNodesEl, EMatrix::Identity(spaceDim, spaceDim));
  }
  //start with face integrals
  for(int iFace = 0; iFace < nFaces; iFace++){
    faceNodes = &(faceNodeMap->at(iFace));
    for(int ip = 0; ip < nIPsFc; ip++){
      shapes = &(fipShapes->at(ip));
      offset = iFace*nIPsFc + ip;
      buffVec = Ds[nIPsEl + offset] * (normals->at(offset)) * dV[nIPsEl + offset];
      for(int iN = 0; iN < nNodesPFc; iN++){
        buffVec2 = buffVec * shapes->at(iN);
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          for(int jN = 0; jN < nNodesPFc; jN++){
            buffVec3 = buffVec2 * shapes->at(jN);
            for(int d = 0; d < spaceDim; d++){
              op(lenU + lenQ + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, lenU + (faceNodes->at(jN)*spaceDim + d)*nDOFsPerNode + ndof) -= buffVec3[d];
              op(faceNodes->at(iN)*nDOFsPerNode + ndof, lenU + (faceNodes->at(jN)*spaceDim + d)*nDOFsPerNode + ndof) -= buffVec3[d];
            }
          }
        }
      }
    }
  }
  //go into bulk integral
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    for(int iN = 0; iN < nNodesEl; iN++){
      buffVec = Ds[ip]*invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), derivShapes->at(iN).size())*dV[ip];
      for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
        for(int jN = 0; jN < nNodesEl; jN++){
          buffVec2 = buffVec * shapes->at(jN);
          for(int d = 0; d < spaceDim; d++){
            op(iN*nDOFsPerNode + ndof, lenU + (jN*spaceDim + d)*nDOFsPerNode + ndof) += buffVec2[d];
          }
        }
      }
    }
  }
};//assemble

}//hfox
